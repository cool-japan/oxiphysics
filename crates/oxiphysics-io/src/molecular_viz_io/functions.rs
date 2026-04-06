//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Covalent radius (Å) and CPK color (r, g, b) for common elements.
pub(super) fn element_data(symbol: &str) -> (f64, [f32; 3]) {
    match symbol.trim().to_uppercase().as_str() {
        "H" => (0.31, [1.00, 1.00, 1.00]),
        "C" => (0.76, [0.50, 0.50, 0.50]),
        "N" => (0.71, [0.13, 0.18, 0.80]),
        "O" => (0.66, [0.80, 0.13, 0.13]),
        "S" => (1.05, [1.00, 0.78, 0.20]),
        "P" => (1.07, [1.00, 0.50, 0.00]),
        "F" => (0.57, [0.56, 0.88, 0.31]),
        "CL" => (1.02, [0.12, 0.94, 0.12]),
        "BR" => (1.20, [0.65, 0.16, 0.16]),
        "I" => (1.39, [0.58, 0.00, 0.58]),
        "NA" => (1.66, [0.67, 0.36, 0.95]),
        "K" => (2.03, [0.56, 0.25, 0.83]),
        "CA" => (1.76, [0.24, 1.00, 0.00]),
        "MG" => (1.41, [0.54, 1.00, 0.00]),
        "ZN" => (1.22, [0.49, 0.50, 0.69]),
        "FE" => (1.32, [0.88, 0.40, 0.20]),
        "CU" => (1.32, [0.78, 0.50, 0.20]),
        _ => (1.50, [0.80, 0.80, 0.80]),
    }
}
/// Van der Waals radius (Å) for common elements.
pub(super) fn vdw_radius(symbol: &str) -> f64 {
    match symbol.trim().to_uppercase().as_str() {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.80,
        "P" => 1.80,
        "F" => 1.47,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        "NA" => 2.27,
        "K" => 2.75,
        "CA" => 2.31,
        "MG" => 1.73,
        "ZN" => 1.39,
        "FE" => 2.00,
        "CU" => 1.40,
        _ => 1.70,
    }
}
/// Sanitize a name for use as a PyMOL or Tcl identifier.
pub(super) fn sanitize_name(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}
/// Escape a string for use inside a JSON string value.
pub(super) fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
}
/// Format an element symbol for use in MDL Molfile (3-char field).
pub(super) fn format_element_symbol(sym: &str) -> String {
    let upper = sym.trim().to_uppercase();
    match upper.as_str() {
        "C" => "C".to_string(),
        "H" => "H".to_string(),
        "N" => "N".to_string(),
        "O" => "O".to_string(),
        "S" => "S".to_string(),
        "P" => "P".to_string(),
        "F" => "F".to_string(),
        "CL" => "Cl".to_string(),
        "BR" => "Br".to_string(),
        "I" => "I".to_string(),
        s => {
            let mut chars = s.chars();
            match chars.next() {
                None => "C".to_string(),
                Some(first) => {
                    let rest: String = chars.collect::<String>().to_lowercase();
                    format!("{}{}", first, rest)
                }
            }
        }
    }
}
/// Return the atomic number for common elements.
pub(super) fn atomic_number(symbol: &str) -> u8 {
    match symbol.trim().to_uppercase().as_str() {
        "H" => 1,
        "HE" => 2,
        "LI" => 3,
        "BE" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "NE" => 10,
        "NA" => 11,
        "MG" => 12,
        "AL" => 13,
        "SI" => 14,
        "P" => 15,
        "S" => 16,
        "CL" => 17,
        "AR" => 18,
        "K" => 19,
        "CA" => 20,
        "FE" => 26,
        "CU" => 29,
        "ZN" => 30,
        "BR" => 35,
        "I" => 53,
        _ => 6,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecular_viz_io::AtomGroup;
    use crate::molecular_viz_io::ChemicalJsonWriter;
    use crate::molecular_viz_io::ColorScheme;
    use crate::molecular_viz_io::MolAtom;
    use crate::molecular_viz_io::MolBond;
    use crate::molecular_viz_io::MolecularScene;
    use crate::molecular_viz_io::MolfileWriter;
    use crate::molecular_viz_io::PymolScriptWriter;
    use crate::molecular_viz_io::VmdRepresentation;
    use crate::molecular_viz_io::VmdScriptWriter;
    use crate::molecular_viz_io::XyzFrame;
    use crate::molecular_viz_io::XyzTrajectoryWriter;
    fn water_scene() -> MolecularScene {
        let mut scene = MolecularScene::new("water");
        let o = scene.add_atom(MolAtom::new(0, "O", [0.000, 0.000, 0.000]));
        let h1 = scene.add_atom(MolAtom::new(1, "H", [0.960, 0.000, 0.000]));
        let h2 = scene.add_atom(MolAtom::new(2, "H", [-0.240, 0.927, 0.000]));
        scene.add_bond(MolBond::single(o, h1));
        scene.add_bond(MolBond::single(o, h2));
        scene
    }
    fn benzene_scene() -> MolecularScene {
        let mut scene = MolecularScene::new("benzene");
        let r = 1.4_f64;
        let positions: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let theta = i as f64 * std::f64::consts::PI / 3.0;
                [r * theta.cos(), r * theta.sin(), 0.0]
            })
            .collect();
        for (i, pos) in positions.iter().enumerate() {
            scene.add_atom(MolAtom::new(i, "C", *pos));
        }
        for i in 0..6 {
            scene.add_bond(MolBond::aromatic(i, (i + 1) % 6));
        }
        scene
    }
    #[test]
    fn test_mol_atom_new_oxygen() {
        let a = MolAtom::new(0, "O", [0.0, 0.0, 0.0]);
        assert_eq!(a.element, "O");
        assert!(a.radius > 0.0);
        assert_eq!(a.index, 0);
    }
    #[test]
    fn test_mol_atom_cpk_color_hydrogen() {
        let h = MolAtom::new(0, "H", [0.0, 0.0, 0.0]);
        let c = h.cpk_color();
        assert!((c[0] - 1.0).abs() < 1e-5);
        assert!((c[1] - 1.0).abs() < 1e-5);
        assert!((c[2] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_mol_atom_cpk_color_oxygen() {
        let o = MolAtom::new(0, "O", [0.0, 0.0, 0.0]);
        let c = o.cpk_color();
        assert!(c[0] > c[1], "O should be reddish");
        assert!(c[0] > c[2], "O should be reddish");
    }
    #[test]
    fn test_mol_atom_covalent_radius_carbon() {
        let c = MolAtom::new(0, "C", [0.0, 0.0, 0.0]);
        assert!((c.covalent_radius() - 0.76).abs() < 0.01);
    }
    #[test]
    fn test_mol_atom_unknown_element_defaults() {
        let a = MolAtom::new(0, "Xe", [0.0, 0.0, 0.0]);
        assert!(a.radius > 0.0);
        assert_eq!(a.element, "Xe");
    }
    #[test]
    fn test_bond_orders() {
        assert_eq!(MolBond::single(0, 1).order.as_int(), 1);
        assert_eq!(MolBond::double(0, 1).order.as_int(), 2);
        assert_eq!(MolBond::triple(0, 1).order.as_int(), 3);
        assert_eq!(MolBond::aromatic(0, 1).order.as_float(), 1.5);
    }
    #[test]
    fn test_bond_aromatic_flag() {
        let b = MolBond::aromatic(0, 1);
        assert!(b.is_aromatic);
        assert!(b.in_ring);
    }
    #[test]
    fn test_bond_single_not_aromatic() {
        let b = MolBond::single(0, 1);
        assert!(!b.is_aromatic);
    }
    #[test]
    fn test_atom_group_add_and_len() {
        let mut grp = AtomGroup::new("test");
        assert!(grp.is_empty());
        grp.add(0);
        grp.add(1);
        assert_eq!(grp.len(), 2);
        assert!(!grp.is_empty());
    }
    #[test]
    fn test_atom_group_from_indices() {
        let grp = AtomGroup::from_indices("selected", vec![1, 3, 5]);
        assert_eq!(grp.len(), 3);
        assert_eq!(grp.name, "selected");
    }
    #[test]
    fn test_scene_add_atoms_and_bonds() {
        let scene = water_scene();
        assert_eq!(scene.natoms(), 3);
        assert_eq!(scene.nbonds(), 2);
    }
    #[test]
    fn test_scene_centroid_water() {
        let scene = water_scene();
        let c = scene.centroid();
        assert!(c[1] > 0.0);
    }
    #[test]
    fn test_scene_bounding_box() {
        let scene = water_scene();
        let (mn, mx) = scene.bounding_box();
        assert!(mn[0] <= 0.0);
        assert!(mx[0] >= 0.0);
        assert!(mn[1] <= 0.0);
        assert!(mx[1] >= 0.0);
    }
    #[test]
    fn test_scene_empty_centroid() {
        let scene = MolecularScene::new("empty");
        let c = scene.centroid();
        assert_eq!(c, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_scene_auto_bonds() {
        let mut scene = MolecularScene::new("water");
        scene.add_atom(MolAtom::new(0, "O", [0.0, 0.0, 0.0]));
        scene.add_atom(MolAtom::new(1, "H", [0.96, 0.0, 0.0]));
        scene.add_atom(MolAtom::new(2, "H", [-0.24, 0.927, 0.0]));
        scene.auto_bonds();
        assert!(scene.nbonds() >= 2);
    }
    #[test]
    fn test_scene_group_by_chain() {
        let mut scene = water_scene();
        scene.group_by_chain();
        assert!(scene.groups.contains_key("chain_A"));
    }
    #[test]
    fn test_scene_set_property() {
        let mut scene = water_scene();
        scene.set_property("energy", -75.123);
        assert!((scene.properties["energy"] - (-75.123)).abs() < 1e-10);
    }
    #[test]
    fn test_scene_set_metadata() {
        let mut scene = water_scene();
        scene.set_metadata("basis", "6-31G*");
        assert_eq!(scene.metadata["basis"], "6-31G*");
    }
    #[test]
    fn test_scene_unit_cell() {
        let mut scene = water_scene();
        scene.set_unit_cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0]);
        assert!(scene.unit_cell.is_some());
    }
    #[test]
    fn test_scene_atom_accessor() {
        let scene = water_scene();
        assert!(scene.atom(0).is_some());
        assert!(scene.atom(99).is_none());
        assert_eq!(scene.atom(0).unwrap().element, "O");
    }
    #[test]
    fn test_pymol_script_contains_show_sticks() {
        let scene = water_scene();
        let writer = PymolScriptWriter::new();
        let script = writer.write_scene(&scene);
        assert!(
            script.contains("show sticks"),
            "Script must contain 'show sticks'"
        );
    }
    #[test]
    fn test_pymol_script_contains_bg_color() {
        let scene = water_scene();
        let writer = PymolScriptWriter::new();
        let script = writer.write_scene(&scene);
        assert!(script.contains("bg_color"));
    }
    #[test]
    fn test_pymol_script_show_spheres_flag() {
        let scene = water_scene();
        let writer = PymolScriptWriter {
            show_spheres: true,
            ..PymolScriptWriter::new()
        };
        let script = writer.write_scene(&scene);
        assert!(script.contains("show spheres"));
    }
    #[test]
    fn test_pymol_script_no_spheres_by_default() {
        let scene = water_scene();
        let writer = PymolScriptWriter::new();
        let script = writer.write_scene(&scene);
        assert!(!script.contains("show spheres"));
    }
    #[test]
    fn test_pymol_script_ray_command() {
        let scene = water_scene();
        let writer = PymolScriptWriter {
            add_ray: true,
            ..PymolScriptWriter::new()
        };
        let script = writer.write_scene(&scene);
        assert!(script.contains("\nray\n") || script.ends_with("ray\n"));
    }
    #[test]
    fn test_pymol_bfactor_spectrum() {
        let mut scene = water_scene();
        for a in &mut scene.atoms {
            a.scalar_property = 10.0;
        }
        let writer = PymolScriptWriter::new();
        let script = writer.write_bfactor_spectrum(&scene);
        assert!(script.contains("spectrum"));
    }
    #[test]
    fn test_pymol_per_atom_colors() {
        let scene = water_scene();
        let writer = PymolScriptWriter::new();
        let script = writer.write_per_atom_colors(&scene);
        assert!(script.contains("set_color atom_color_0"));
    }
    #[test]
    fn test_pymol_color_scheme_by_chain() {
        let scene = water_scene();
        let writer = PymolScriptWriter {
            color_scheme: ColorScheme::ByChain,
            ..PymolScriptWriter::new()
        };
        let script = writer.write_scene(&scene);
        assert!(script.contains("util.cbc"));
    }
    #[test]
    fn test_pymol_group_selection_written() {
        let mut scene = water_scene();
        let mut grp = AtomGroup::from_indices("heavy", vec![0]);
        grp.color = Some([0.8, 0.2, 0.1]);
        scene.add_group(grp);
        let writer = PymolScriptWriter::new();
        let script = writer.write_scene(&scene);
        assert!(script.contains("heavy"));
    }
    #[test]
    fn test_vmd_script_mol_load() {
        let scene = water_scene();
        let writer = VmdScriptWriter::new();
        let script = writer.write_scene(&scene);
        assert!(script.contains("mol load xyz"));
    }
    #[test]
    fn test_vmd_script_addrep() {
        let scene = water_scene();
        let writer = VmdScriptWriter::new();
        let script = writer.write_scene(&scene);
        assert!(script.contains("mol addrep"));
    }
    #[test]
    fn test_vmd_script_representation_licorice() {
        let scene = water_scene();
        let writer = VmdScriptWriter {
            representation: VmdRepresentation::Licorice,
            ..VmdScriptWriter::new()
        };
        let script = writer.write_scene(&scene);
        assert!(script.contains("Licorice"));
    }
    #[test]
    fn test_vmd_script_representation_cpk() {
        let scene = water_scene();
        let writer = VmdScriptWriter {
            representation: VmdRepresentation::Cpk,
            ..VmdScriptWriter::new()
        };
        let script = writer.write_scene(&scene);
        assert!(script.contains("CPK"));
    }
    #[test]
    fn test_vmd_trajectory_script() {
        let scene = water_scene();
        let frame1: Vec<MolAtom> = scene.atoms.clone();
        let frame2: Vec<MolAtom> = scene
            .atoms
            .iter()
            .map(|a| {
                let mut b = a.clone();
                b.position[0] += 0.01;
                b
            })
            .collect();
        let writer = VmdScriptWriter::new();
        let script = writer.write_trajectory_script(&[frame1, frame2], "water_traj");
        assert!(script.contains("mol addfile"));
    }
    #[test]
    fn test_vmd_backbone_representation() {
        let scene = water_scene();
        let writer = VmdScriptWriter {
            show_backbone: true,
            ..VmdScriptWriter::new()
        };
        let script = writer.write_scene(&scene);
        assert!(script.contains("NewRibbons"));
    }
    #[test]
    fn test_chemical_json_contains_atoms() {
        let scene = water_scene();
        let writer = ChemicalJsonWriter::new();
        let json = writer.write_scene(&scene);
        assert!(json.contains("\"atoms\""));
    }
    #[test]
    fn test_chemical_json_contains_bonds() {
        let scene = water_scene();
        let writer = ChemicalJsonWriter::new();
        let json = writer.write_scene(&scene);
        assert!(json.contains("\"bonds\""));
    }
    #[test]
    fn test_chemical_json_contains_elements() {
        let scene = water_scene();
        let writer = ChemicalJsonWriter::new();
        let json = writer.write_scene(&scene);
        assert!(json.contains("\"elements\""));
    }
    #[test]
    fn test_chemical_json_oxygen_atomic_number() {
        assert_eq!(atomic_number("O"), 8);
        assert_eq!(atomic_number("C"), 6);
        assert_eq!(atomic_number("H"), 1);
        assert_eq!(atomic_number("N"), 7);
    }
    #[test]
    fn test_chemical_json_coords_3d() {
        let scene = water_scene();
        let writer = ChemicalJsonWriter::new();
        let json = writer.write_scene(&scene);
        assert!(json.contains("\"3d\""));
    }
    #[test]
    fn test_chemical_json_benzene_bond_orders() {
        let scene = benzene_scene();
        let writer = ChemicalJsonWriter::new();
        let json = writer.write_scene(&scene);
        assert!(json.contains("\"order\""));
    }
    #[test]
    fn test_chemical_json_properties() {
        let mut scene = water_scene();
        scene.set_property("totalEnergy", -75.1234);
        let writer = ChemicalJsonWriter::new();
        let json = writer.write_scene(&scene);
        assert!(json.contains("\"totalEnergy\""));
    }
    #[test]
    fn test_molfile_v2000_header() {
        let scene = water_scene();
        let writer = MolfileWriter::new();
        let mol = writer.write_scene(&scene);
        assert!(mol.contains("V2000"));
    }
    #[test]
    fn test_molfile_m_end() {
        let scene = water_scene();
        let writer = MolfileWriter::new();
        let mol = writer.write_scene(&scene);
        assert!(mol.contains("M  END"));
    }
    #[test]
    fn test_molfile_atom_count_in_header() {
        let scene = water_scene();
        let writer = MolfileWriter::new();
        let mol = writer.write_scene(&scene);
        assert!(mol.contains("  3  2"));
    }
    #[test]
    fn test_molfile_sdf_mode_separator() {
        let scene = water_scene();
        let writer = MolfileWriter {
            sdf_mode: true,
            ..MolfileWriter::new()
        };
        let mol = writer.write_scene(&scene);
        assert!(mol.contains("$$$$"));
    }
    #[test]
    fn test_molfile_sdf_multiple_scenes() {
        let scenes = vec![water_scene(), benzene_scene()];
        let writer = MolfileWriter::new();
        let sdf = writer.write_sdf(&scenes);
        assert_eq!(sdf.matches("$$$$").count(), 2);
    }
    #[test]
    fn test_molfile_element_symbols() {
        assert_eq!(format_element_symbol("C"), "C");
        assert_eq!(format_element_symbol("CL"), "Cl");
        assert_eq!(format_element_symbol("BR"), "Br");
        assert_eq!(format_element_symbol("H"), "H");
    }
    #[test]
    fn test_xyz_frame_natoms() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let elements = vec!["H".to_string(), "H".to_string()];
        let frame = XyzFrame::new(positions, elements);
        assert_eq!(frame.natoms(), 2);
    }
    #[test]
    fn test_xyz_writer_frame_header() {
        let scene = water_scene();
        let frame = XyzFrame::from_scene(&scene, 0, 0.0);
        let writer = XyzTrajectoryWriter::new();
        let out = writer.write_frame(&frame);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0].trim(), "3");
    }
    #[test]
    fn test_xyz_writer_step_in_comment() {
        let scene = water_scene();
        let frame = XyzFrame::from_scene(&scene, 42, 1.23);
        let writer = XyzTrajectoryWriter::new();
        let out = writer.write_frame(&frame);
        let comment = out.lines().nth(1).unwrap();
        assert!(comment.contains("step=42"));
        assert!(comment.contains("time=1.23"));
    }
    #[test]
    fn test_xyz_writer_energy_in_comment() {
        let scene = water_scene();
        let mut frame = XyzFrame::from_scene(&scene, 0, 0.0);
        frame.total_energy = Some(-75.1234);
        let writer = XyzTrajectoryWriter::new();
        let out = writer.write_frame(&frame);
        assert!(out.contains("energy="));
    }
    #[test]
    fn test_xyz_writer_element_symbols() {
        let scene = water_scene();
        let frame = XyzFrame::from_scene(&scene, 0, 0.0);
        let writer = XyzTrajectoryWriter::new();
        let out = writer.write_frame(&frame);
        assert!(out.contains("O  ") || out.contains("O "));
        assert!(out.contains("H  ") || out.contains("H "));
    }
    #[test]
    fn test_xyz_writer_multi_frame() {
        let scene = water_scene();
        let frames: Vec<XyzFrame> = (0..3)
            .map(|i| XyzFrame::from_scene(&scene, i, i as f64 * 0.01))
            .collect();
        let writer = XyzTrajectoryWriter::new();
        let out = writer.write_trajectory(&frames);
        assert_eq!(out.matches("step=").count(), 3);
    }
    #[test]
    fn test_xyz_writer_velocities() {
        let mut frame = XyzFrame::new(vec![[0.0, 0.0, 0.0]], vec!["C".to_string()]);
        frame.velocities = Some(vec![[1.0, 2.0, 3.0]]);
        let writer = XyzTrajectoryWriter {
            include_velocities: true,
            ..XyzTrajectoryWriter::new()
        };
        let out = writer.write_frame(&frame);
        let atom_line = out.lines().nth(2).unwrap();
        let parts: Vec<&str> = atom_line.split_whitespace().collect();
        assert!(
            parts.len() >= 7,
            "expected 7+ fields (element + 6 floats), got {}",
            parts.len()
        );
    }
    #[test]
    fn test_xyz_writer_append_frame() {
        let scene = water_scene();
        let frame = XyzFrame::from_scene(&scene, 0, 0.0);
        let writer = XyzTrajectoryWriter::new();
        let mut buf = String::new();
        writer.append_frame(&mut buf, &frame);
        writer.append_frame(&mut buf, &frame);
        assert_eq!(buf.matches("step=").count(), 2);
    }
    #[test]
    fn test_xyz_with_dynamics_writer() {
        let writer = XyzTrajectoryWriter::with_dynamics();
        assert!(writer.include_velocities);
        assert!(writer.include_forces);
    }
    #[test]
    fn test_sanitize_name() {
        assert_eq!(sanitize_name("my scene"), "my_scene");
        assert_eq!(sanitize_name("abc-123"), "abc_123");
        assert_eq!(sanitize_name("valid_name"), "valid_name");
    }
    #[test]
    fn test_escape_json() {
        assert_eq!(escape_json("hello"), "hello");
        assert_eq!(escape_json("a\"b"), "a\\\"b");
        assert_eq!(escape_json("a\\b"), "a\\\\b");
    }
    #[test]
    fn test_vdw_radius_values() {
        assert!((vdw_radius("H") - 1.20).abs() < 0.01);
        assert!((vdw_radius("C") - 1.70).abs() < 0.01);
        assert!((vdw_radius("O") - 1.52).abs() < 0.01);
    }
    #[test]
    fn test_element_data_sulfur() {
        let (r, _c) = element_data("S");
        assert!((r - 1.05).abs() < 0.01);
    }
}
