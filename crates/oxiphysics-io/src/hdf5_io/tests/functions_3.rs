//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod microstructure_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_phase_field_2d() {
        let mut f = Hdf5File::new("t.h5");
        let phi = vec![0.0, 0.5, 1.0, 0.5];
        write_phase_field_2d(&mut f, "pf", 0, 2, 2, &phi).unwrap();
        let g = f.open_group("pf").unwrap();
        assert!(g.groups.contains_key("step_0000000000"));
    }
    #[test]
    fn test_write_chemical_potential_field() {
        let mut f = Hdf5File::new("t.h5");
        let mu = vec![0.1, -0.1, 0.2, -0.2];
        write_chemical_potential_field(&mut f, "ch", 2, 2, &mu).unwrap();
        assert_eq!(
            f.open_dataset("ch", "chemical_potential")
                .unwrap()
                .read_f64()
                .unwrap(),
            mu
        );
    }
    #[test]
    fn test_write_percolation_stats() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![0.5, 0.6, 0.7];
        let cs = vec![3.0, 10.0, 100.0];
        let perc = vec![false, false, true];
        write_percolation_stats(&mut f, "perc", &p, &cs, &perc).unwrap();
        let pf = f
            .open_dataset("perc", "percolates")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(pf, vec![0.0, 0.0, 1.0]);
    }
    #[test]
    fn test_write_heat_field() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![300.0, 350.0, 400.0, 350.0];
        write_heat_field(&mut f, "heat", 5, 2, 2, &t).unwrap();
        let g = f.open_group("heat").unwrap();
        assert!(g.groups.contains_key("step_0000000005"));
    }
    #[test]
    fn test_write_concentration_field() {
        let mut f = Hdf5File::new("t.h5");
        let c = vec![0.1, 0.2, 0.3, 0.4];
        write_concentration_field(&mut f, "rd", "A", 0, 2, 2, &c).unwrap();
        let back = f
            .open_dataset("rd/A/step_0000000000", "concentration")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(back, c);
    }
    #[test]
    fn test_write_nucleation_data() {
        let mut f = Hdf5File::new("t.h5");
        let ns = vec![1usize, 5, 10, 20];
        let fg = vec![0.5, -0.2, -1.0, -2.5];
        write_nucleation_data(&mut f, "nucl", &ns, &fg).unwrap();
        assert_eq!(
            f.open_dataset("nucl", "free_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            fg
        );
    }
}
#[cfg(test)]
mod extra_roundtrip_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_soc_matrix_write() {
        let mut f = Hdf5File::new("t.h5");
        let mat = vec![0.0_f64; 32];
        write_soc_matrix(&mut f, "soc", 2, &mat).unwrap();
        assert_eq!(
            f.open_dataset("soc", "soc_matrix")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            32
        );
    }
    #[test]
    fn test_write_gw_corrections() {
        let mut f = Hdf5File::new("t.h5");
        let dft = vec![-3.0, -1.0, 2.0];
        let gw = vec![-3.5, -1.2, 2.3];
        write_gw_corrections(&mut f, "gw2", &dft, &gw).unwrap();
        assert_eq!(
            f.open_dataset("gw2", "gw_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            gw
        );
    }
    #[test]
    fn test_write_force_constants() {
        let mut f = Hdf5File::new("t.h5");
        let ifc = vec![0.1_f64; 36];
        write_force_constants(&mut f, "ifc", 2, &ifc).unwrap();
        assert_eq!(
            f.open_dataset("ifc", "force_constants")
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
        write_dipole_trajectory(&mut f, "dip", &d).unwrap();
        assert_eq!(
            f.open_dataset("dip", "dipole_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }
    #[test]
    fn test_write_magnetic_moments() {
        let mut f = Hdf5File::new("t.h5");
        let m = vec![[0.0, 0.0, 2.0]];
        write_magnetic_moments(&mut f, "mag", &m).unwrap();
        assert_eq!(
            f.open_dataset("mag", "magnetic_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            3
        );
    }
    #[test]
    fn test_write_electron_density() {
        let mut f = Hdf5File::new("t.h5");
        let rho = vec![0.5_f64; 27];
        write_electron_density(&mut f, "dft2", 3, 3, 3, &rho).unwrap();
        assert_eq!(
            f.open_dataset("dft2", "rho_e")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            27
        );
    }
    #[test]
    fn test_write_convex_hull_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let x = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let e = vec![0.0, -0.05, -0.10, -0.05, 0.0];
        let s = vec![true, true, true, true, true];
        write_convex_hull(&mut f, "hull2", &x, &e, &s).unwrap();
        assert_eq!(
            f.open_dataset("hull2", "hull_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            e
        );
    }
    #[test]
    fn test_write_thermochemistry_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![298.15, 500.0, 1000.0];
        let h = vec![-200.0, -195.0, -185.0];
        let s = vec![0.1, 0.15, 0.25];
        let g: Vec<f64> = h
            .iter()
            .zip(t.iter())
            .zip(s.iter())
            .map(|((&hi, &ti), &si)| hi - ti * si)
            .collect();
        write_thermochemistry(&mut f, "tc", &t, &h, &s, &g).unwrap();
        assert_eq!(
            f.open_dataset("tc", "temperatures")
                .unwrap()
                .read_f64()
                .unwrap(),
            t
        );
    }
    #[test]
    fn test_write_optical_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![1.0, 2.0, 3.0, 4.0];
        let eps2 = vec![0.0, 0.1, 0.5, 0.2];
        write_optical_spectrum(&mut f, "opt2", &e, &eps2).unwrap();
        assert_eq!(
            f.open_dataset("opt2", "epsilon2")
                .unwrap()
                .read_f64()
                .unwrap(),
            eps2
        );
    }
    #[test]
    fn test_write_lattice_parameters_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![100.0, 300.0, 500.0];
        let a = vec![5.0, 5.02, 5.05];
        write_lattice_parameters(&mut f, "th", &t, &a, &a, &a).unwrap();
        assert_eq!(f.open_dataset("th", "a").unwrap().read_f64().unwrap(), a);
    }
    #[test]
    fn test_hdf5_file_root_is_slash() {
        let f = Hdf5File::new("t.h5");
        assert_eq!(f.root.name, "/");
    }
    #[test]
    fn test_write_neb_path_barrier() {
        let mut f = Hdf5File::new("t.h5");
        let im = vec![
            vec![[0.0, 0.0, 0.0]],
            vec![[0.5, 0.0, 0.0]],
            vec![[1.0, 0.0, 0.0]],
        ];
        let e = vec![-10.0, -7.0, -10.0];
        write_neb_path(&mut f, "neb2", &im, &e).unwrap();
        let eb = read_neb_energies(&f, "neb2").unwrap();
        assert!((eb[1] - (-7.0)).abs() < 1e-12);
    }
    #[test]
    fn test_write_wannier_centres_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let c = vec![[0.5, 0.5, 0.5]];
        let sp = vec![1.2];
        write_wannier_centres(&mut f, "wan", &c, &sp).unwrap();
        assert_eq!(
            f.open_dataset("wan", "wannier_spreads")
                .unwrap()
                .read_f64()
                .unwrap(),
            sp
        );
    }
}
#[cfg(test)]
mod integration_final_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_convergence_log() {
        let mut f = Hdf5File::new("t.h5");
        write_convergence_log(&mut f, "scf", &[1, 2, 3], &[1e-2, 1e-4, 1e-8], true).unwrap();
        let g = f.open_group("scf").unwrap();
        assert_eq!(g.attributes["converged"], AttrValue::Int32(1));
    }
    #[test]
    fn test_write_geopt_step() {
        let mut f = Hdf5File::new("t.h5");
        let pos = vec![[0.0, 0.0, 0.0]];
        let frc = vec![[0.05, 0.0, 0.0]];
        write_geopt_step(&mut f, "geopt", 0, &pos, &frc, -10.5, 0.05).unwrap();
        let g = f.open_group("geopt").unwrap();
        let sg = g.groups.get("geopt_000000").unwrap();
        assert_eq!(sg.attributes["max_force"], AttrValue::Float64(0.05));
    }
    #[test]
    fn test_write_mulliken_charges() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.3, -0.3, 0.3, -0.3];
        write_mulliken_charges(&mut f, "pop", &q, None).unwrap();
        assert_eq!(
            f.open_dataset("pop", "mulliken_charges")
                .unwrap()
                .read_f64()
                .unwrap(),
            q
        );
    }
    #[test]
    fn test_write_mulliken_with_spin() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.1, -0.1];
        let sd = vec![0.5, -0.5];
        write_mulliken_charges(&mut f, "rad", &q, Some(&sd)).unwrap();
        assert_eq!(
            f.open_dataset("rad", "spin_densities")
                .unwrap()
                .read_f64()
                .unwrap(),
            sd
        );
    }
    #[test]
    fn test_write_parallel_coordinates() {
        let mut f = Hdf5File::new("t.h5");
        let vn = vec!["x".to_string(), "y".to_string(), "z".to_string()];
        let data: Vec<f64> = (0..9).map(|i| i as f64).collect();
        write_parallel_coordinates(&mut f, "pcoord", &vn, &data, 3).unwrap();
        assert_eq!(
            f.open_dataset("pcoord", "data")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            9
        );
    }
    #[test]
    fn test_write_percolation_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![0.4, 0.5, 0.6, 0.7];
        let cs = vec![2.0, 4.0, 15.0, 100.0];
        let perc = vec![false, false, false, true];
        write_percolation_stats(&mut f, "p2", &p, &cs, &perc).unwrap();
        let pf = f
            .open_dataset("p2", "percolates")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(pf[3], 1.0);
        assert_eq!(pf[0], 0.0);
    }
    #[test]
    fn test_write_eos_fit_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![60.0, 65.0, 70.0];
        let e = vec![-4.0, -4.5, -4.3];
        write_eos_fit(&mut f, "eos2", &v, &e, 66.0, -4.5, 150.0, 4.2).unwrap();
        let g = f.open_group("eos2").unwrap();
        assert_eq!(g.attributes["E0"], AttrValue::Float64(-4.5));
    }
    #[test]
    fn test_total_dataset_count_with_groups() {
        let mut f = Hdf5File::new("t.h5");
        write_masses(&mut f, "atoms", &[1.0, 2.0, 3.0]).unwrap();
        write_charges(&mut f, "atoms", &[0.1, -0.1, 0.0]).unwrap();
        write_rdf(&mut f, "obs", &[0.1, 0.2], &[0.5, 1.0]).unwrap();
        let n = total_dataset_count(&f);
        assert!(n >= 4, "expected >=4, got {n}");
    }
    #[test]
    fn test_list_top_level_groups_non_empty() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 1, 0, 0, "OxiPhysics").unwrap();
        write_masses(&mut f, "atoms", &[1.0]).unwrap();
        let groups = list_top_level_groups(&f);
        assert!(groups.contains(&"atoms".to_string()));
        assert!(groups.contains(&"__version__".to_string()));
    }
    #[test]
    fn test_write_read_band_structure_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let k: Vec<f64> = (0..9).map(|i| i as f64 * 0.1).collect();
        let e: Vec<f64> = (0..6).map(|i| i as f64 - 3.0).collect();
        write_band_structure(&mut f, "bs2", 3, 2, &k, &e).unwrap();
        let bk = f
            .open_dataset("bs2", "kpoints")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(bk.len(), 9);
        let be = f
            .open_dataset("bs2", "eigenvalues")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(be, e);
    }
    #[test]
    fn test_write_read_velocities_single() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![[1.0, 2.0, 3.0]];
        write_velocities(&mut f, "single_v", &v).unwrap();
        let b = read_velocities(&f, "single_v").unwrap();
        assert_eq!(b, v);
    }
    #[test]
    fn test_write_read_box_vectors_ortho() {
        let mut f = Hdf5File::new("t.h5");
        let mut bv = [0.0_f64; 9];
        bv[0] = 5.0;
        bv[4] = 5.0;
        bv[8] = 5.0;
        write_box_vectors(&mut f, "box2", &bv).unwrap();
        let b = read_box_vectors(&f, "box2").unwrap();
        assert!((b[0] - 5.0).abs() < 1e-12);
        assert!((b[4] - 5.0).abs() < 1e-12);
    }
}
#[cfg(test)]
mod graph_frag_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_molecular_graph() {
        let mut f = Hdf5File::new("t.h5");
        let bonds = vec![(0usize, 1, 1.0), (1, 2, 2.0), (2, 3, 1.5)];
        write_molecular_graph(&mut f, "graph", &bonds).unwrap();
        assert_eq!(
            f.open_dataset("graph", "bond_order")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![1.0, 2.0, 1.5]
        );
    }
    #[test]
    fn test_write_fragmentation_map() {
        let mut f = Hdf5File::new("t.h5");
        let af = vec![0usize, 0, 1, 1, 2];
        write_fragmentation_map(&mut f, "frag", &af).unwrap();
        let b = f
            .open_dataset("frag", "atom_fragment")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(b, vec![0.0, 0.0, 1.0, 1.0, 2.0]);
    }
    #[test]
    fn test_write_fragment_masses() {
        let mut f = Hdf5File::new("t.h5");
        write_fragment_masses(&mut f, "frm", &[18.015, 44.010]).unwrap();
        assert_eq!(
            f.open_dataset("frm", "fragment_masses")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![18.015, 44.010]
        );
    }
    #[test]
    fn test_write_read_md_metadata_full() {
        let mut f = Hdf5File::new("t.h5");
        let m = MdRunMetadata {
            title: "NaCl Melt".into(),
            n_atoms: 512,
            n_steps: 500_000,
            dt_fs: 1.0,
            temperature_k: 1200.0,
            pressure_bar: 1.0,
            ensemble: "NVT".into(),
            force_field: "Born-Huggins-Meyer".into(),
        };
        m.write_to(&mut f).unwrap();
        let b = MdRunMetadata::read_from(&f).unwrap();
        assert_eq!(b.force_field, "Born-Huggins-Meyer");
        assert!((b.temperature_k - 1200.0).abs() < 1e-12);
    }
    #[test]
    fn test_write_solvation_energies_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let solvents = vec!["water".to_string(), "dmso".to_string()];
        let dg = vec![-4.5, -3.0];
        write_solvation_energies(&mut f, "solv2", &solvents, &dg).unwrap();
        assert_eq!(
            f.open_dataset("solv2", "dG_solv_kcal_mol")
                .unwrap()
                .read_f64()
                .unwrap(),
            dg
        );
    }
    #[test]
    fn test_write_fukui_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let fp = vec![0.3, 0.4, 0.3];
        let fm = vec![0.2, 0.5, 0.3];
        let fz = vec![0.25, 0.45, 0.30];
        write_fukui_functions(&mut f, "fuk2", &fp, &fm, &fz).unwrap();
        assert_eq!(
            f.open_dataset("fuk2", "f_minus")
                .unwrap()
                .read_f64()
                .unwrap(),
            fm
        );
    }
    #[test]
    fn test_write_convergence_log_not_converged() {
        let mut f = Hdf5File::new("t.h5");
        write_convergence_log(
            &mut f,
            "scf2",
            &[1, 2, 3, 4, 5],
            &[1e-1, 1e-2, 1e-3, 1e-4, 1e-5],
            false,
        )
        .unwrap();
        let g = f.open_group("scf2").unwrap();
        assert_eq!(g.attributes["converged"], AttrValue::Int32(0));
    }
    #[test]
    fn test_write_adsorption_energies_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        write_adsorption_energies(&mut f, "ads2", &["bridge", "hollow"], &[-1.8, -2.1]).unwrap();
        let g = f.open_group("ads2").unwrap();
        assert_eq!(g.attributes["hollow"], AttrValue::Float64(-2.1));
    }
    #[test]
    fn test_write_rate_constants() {
        let mut f = Hdf5File::new("t.h5");
        let rxn = vec!["R1".to_string(), "R2".to_string()];
        let kf = vec![1e6, 5e4];
        let kr = vec![1e2, 1e3];
        write_rate_constants(&mut f, "micro", &rxn, &kf, &kr).unwrap();
        assert_eq!(
            f.open_dataset("micro", "k_forward")
                .unwrap()
                .read_f64()
                .unwrap(),
            kf
        );
    }
    #[test]
    fn test_write_coverage_profile_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 1.0, 2.0, 3.0];
        let c = vec![0.0, 0.3, 0.7, 0.9];
        write_coverage_profile(&mut f, "kin", &t, &c, "CO2").unwrap();
        assert_eq!(
            f.open_dataset("kin/CO2", "coverage")
                .unwrap()
                .read_f64()
                .unwrap(),
            c
        );
    }
}
#[cfg(test)]
mod transport_bulk_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_ionic_conductivity() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![300.0, 400.0, 500.0];
        let s = vec![0.01, 0.1, 1.0];
        write_ionic_conductivity(&mut f, "cond", &t, &s).unwrap();
        assert_eq!(
            f.open_dataset("cond", "conductivity_S_m")
                .unwrap()
                .read_f64()
                .unwrap(),
            s
        );
    }
    #[test]
    fn test_write_heat_capacity() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![200.0, 298.0, 400.0];
        let cp = vec![20.0, 25.0, 28.0];
        write_heat_capacity(&mut f, "thcp", &t, &cp).unwrap();
        assert_eq!(
            f.open_dataset("thcp", "Cp_J_mol_K")
                .unwrap()
                .read_f64()
                .unwrap(),
            cp
        );
    }
    #[test]
    fn test_write_cluster_expansion_eci_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let sz = vec![1usize, 2, 3, 4];
        let eci = vec![-0.5, 0.1, 0.05, 0.01];
        write_cluster_expansion_eci(&mut f, "ce2", &sz, &eci).unwrap();
        assert_eq!(
            f.open_dataset("ce2", "eci_values")
                .unwrap()
                .read_f64()
                .unwrap(),
            eci
        );
    }
    #[test]
    fn test_write_nucleation_data_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let ns = vec![1usize, 5, 10];
        let fg = vec![2.0, -0.5, -3.0];
        write_nucleation_data(&mut f, "nucl2", &ns, &fg).unwrap();
        assert_eq!(
            f.open_dataset("nucl2", "free_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            fg
        );
    }
    #[test]
    fn test_write_enumerated_structures() {
        let mut f = Hdf5File::new("t.h5");
        let ids = vec![0usize, 1, 2, 3];
        let x = vec![0.0, 0.25, 0.75, 1.0];
        let ef = vec![0.0, -0.05, -0.05, 0.0];
        write_enumerated_structures(&mut f, "enum2", &ids, &x, &ef).unwrap();
        assert_eq!(
            f.open_dataset("enum2", "concentrations")
                .unwrap()
                .read_f64()
                .unwrap(),
            x
        );
    }
}
#[cfg(test)]
mod scalar_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_read_scalar_with_units() {
        let mut f = Hdf5File::new("t.h5");
        write_scalar_with_units(&mut f, "props", "band_gap", 1.12, "eV").unwrap();
        let v = read_scalar(&f, "props", "band_gap").unwrap();
        assert!((v - 1.12).abs() < 1e-12);
    }
    #[test]
    fn test_read_scalar_error_on_empty() {
        let mut f = Hdf5File::new("t.h5");
        f.create_group("empty").unwrap();
        let _ = f.create_dataset("empty", "val", vec![0], Hdf5Dtype::Float64);
        let result = read_scalar(&f, "empty", "val");
        assert!(result.is_err());
    }
    #[test]
    fn test_write_ionic_conductivity_single_point() {
        let mut f = Hdf5File::new("t.h5");
        write_ionic_conductivity(&mut f, "c1", &[298.15], &[1e-4]).unwrap();
        let s = f
            .open_dataset("c1", "conductivity_S_m")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(s, vec![1e-4]);
    }
    #[test]
    fn test_hdf5_file_is_not_locked_by_default() {
        let f = Hdf5File::new("t.h5");
        assert!(!f.is_locked());
    }
}
#[cfg(test)]
mod root_helper_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_root_group_count_empty() {
        let f = Hdf5File::new("t.h5");
        assert_eq!(root_group_count(&f), 0);
    }
    #[test]
    fn test_is_empty_file() {
        let f = Hdf5File::new("t.h5");
        assert!(is_empty_file(&f));
    }
}
