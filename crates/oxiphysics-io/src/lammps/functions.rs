//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::{Error, Result};

#[cfg(test)]
mod tests {

    use crate::LammpsAtom;
    use crate::LammpsDumpReader;
    use crate::LammpsDumpWriter;

    use crate::lammps::LammpsBox;
    use crate::lammps::LammpsCompute;
    use crate::lammps::LammpsDataReader;

    use crate::lammps::LammpsDataWriter;
    use crate::lammps::LammpsFix;
    use crate::lammps::LammpsHarmonicAngle;
    use crate::lammps::LammpsHarmonicBond;
    use crate::lammps::LammpsInputScript;
    use crate::lammps::LammpsLJPair;
    use crate::lammps::LammpsLJTable;
    use crate::lammps::LammpsLogParser;
    use crate::lammps::LammpsPairStyle;

    use oxiphysics_core::math::Vec3;
    #[test]
    fn test_lammps_dump_write() {
        let atoms = vec![
            LammpsAtom {
                id: 1,
                type_id: 1,
                position: Vec3::new(0.0, 0.0, 0.0),
                velocity: Vec3::new(1.0, 0.0, 0.0),
            },
            LammpsAtom {
                id: 2,
                type_id: 2,
                position: Vec3::new(1.0, 1.0, 1.0),
                velocity: Vec3::new(0.0, 1.0, 0.0),
            },
        ];
        let bounds = [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0]];
        let mut buf = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, 100, bounds, &atoms).unwrap();
        let content = String::from_utf8(buf).unwrap();
        assert!(content.contains("ITEM: TIMESTEP"));
        assert!(content.contains("100"));
        assert!(content.contains("ITEM: NUMBER OF ATOMS"));
        assert!(content.contains("2"));
        assert!(content.contains("ITEM: ATOMS id type x y z vx vy vz"));
    }
    #[test]
    fn test_lammps_dump_write_read() {
        let atoms = vec![
            LammpsAtom {
                id: 1,
                type_id: 1,
                position: Vec3::new(0.1, 0.2, 0.3),
                velocity: Vec3::new(0.0, 0.0, 0.0),
            },
            LammpsAtom {
                id: 2,
                type_id: 1,
                position: Vec3::new(1.0, 2.0, 3.0),
                velocity: Vec3::new(0.1, 0.2, 0.3),
            },
            LammpsAtom {
                id: 3,
                type_id: 2,
                position: Vec3::new(4.0, 5.0, 6.0),
                velocity: Vec3::new(-0.1, -0.2, -0.3),
            },
            LammpsAtom {
                id: 4,
                type_id: 2,
                position: Vec3::new(7.5, 8.5, 9.5),
                velocity: Vec3::new(1.0, 0.0, 0.0),
            },
            LammpsAtom {
                id: 5,
                type_id: 3,
                position: Vec3::new(2.2, 3.3, 4.4),
                velocity: Vec3::new(0.0, 1.0, 0.0),
            },
        ];
        let bounds = [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0]];
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, 42, bounds, &atoms).unwrap();
        let cursor = std::io::Cursor::new(buf);
        let (timestep, read_atoms) = LammpsDumpReader::read_frame(cursor).unwrap();
        assert_eq!(timestep, 42);
        assert_eq!(read_atoms.len(), 5);
        assert_eq!(read_atoms[0].id, 1);
        assert_eq!(read_atoms[1].id, 2);
        assert_eq!(read_atoms[2].id, 3);
        assert_eq!(read_atoms[3].id, 4);
        assert_eq!(read_atoms[4].id, 5);
        assert!((read_atoms[0].position.x - 0.1).abs() < 1e-10);
        assert!((read_atoms[0].position.y - 0.2).abs() < 1e-10);
        assert!((read_atoms[0].position.z - 0.3).abs() < 1e-10);
        assert!((read_atoms[1].position.x - 1.0).abs() < 1e-10);
        assert!((read_atoms[2].position.y - 5.0).abs() < 1e-10);
        assert!((read_atoms[3].position.z - 9.5).abs() < 1e-10);
        assert!((read_atoms[4].position.x - 2.2).abs() < 1e-10);
    }
    #[test]
    fn test_data_writer_header_contains_atoms_count() {
        let h = LammpsDataWriter::write_header(10, 0, [0.0; 3], [5.0, 5.0, 5.0]);
        assert!(h.contains("10 atoms"), "header: {h}");
        assert!(h.contains("LAMMPS data file"));
    }
    #[test]
    fn test_data_writer_box_bounds_in_header() {
        let h = LammpsDataWriter::write_header(1, 0, [-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]);
        assert!(
            h.contains("-1 1 xlo xhi") || h.contains("-1 1"),
            "header should contain box bounds: {h}"
        );
    }
    #[test]
    fn test_data_writer_atomic_section() {
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let masses = vec![12.0, 16.0];
        let types = vec![1u32, 2u32];
        let s = LammpsDataWriter::write_atoms_atomic(&pos, &masses, &types);
        assert!(s.contains("Atoms"), "section header missing: {s}");
        assert!(
            s.contains("1 1 0 0 0") || s.contains("1 1 0.0"),
            "atom 1 missing: {s}"
        );
        assert!(
            s.contains("2 2 1 2 3") || s.contains("2 2 1.0"),
            "atom 2 missing: {s}"
        );
    }
    #[test]
    fn test_data_writer_bonds_section() {
        let bonds = vec![(1u32, 1u32, 2u32), (1, 2, 3)];
        let s = LammpsDataWriter::write_bonds(&bonds);
        assert!(s.contains("Bonds"));
        assert!(s.contains("1 1 1 2") || s.contains("1 1 2"));
    }
    #[test]
    fn test_data_roundtrip() {
        let pos = vec![[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let masses = vec![12.0_f64, 16.0];
        let types = vec![1u32, 2u32];
        let s = LammpsDataWriter::write_complete(&pos, &masses, &types, [0.0; 3], [10.0; 3]);
        let reader = LammpsDataReader::from_str(&s).unwrap();
        assert_eq!(reader.positions().len(), 2);
        assert!((reader.positions()[0][0] - 1.0).abs() < 1e-6);
        assert!((reader.positions()[1][2] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_data_reader_box_bounds() {
        let s = LammpsDataWriter::write_complete(
            &[[0.5, 0.5, 0.5]],
            &[12.0],
            &[1],
            [-5.0, -5.0, -5.0],
            [5.0, 5.0, 5.0],
        );
        let reader = LammpsDataReader::from_str(&s).unwrap();
        let (lo, hi) = reader.box_bounds();
        assert!((lo[0] - (-5.0)).abs() < 1e-6, "lo[0]={}", lo[0]);
        assert!((hi[0] - 5.0).abs() < 1e-6, "hi[0]={}", hi[0]);
    }
    #[test]
    fn test_minimize_script_contains_minimize() {
        let s = LammpsInputScript::minimize_script(1e-4, 1e-8);
        assert!(s.contains("minimize"), "script: {s}");
        assert!(s.contains("0.0001") || s.contains("1e-4") || s.contains("0.00000001"));
    }
    #[test]
    fn test_npt_script_contains_npt() {
        let s = LammpsInputScript::npt_script(300.0, 1.0, 100.0, 1000.0, 10000);
        assert!(s.contains("npt"), "script: {s}");
        assert!(s.contains("300"));
        assert!(s.contains("10000"));
    }
    #[test]
    fn test_write_atoms_full() {
        let pos = vec![[1.0, 2.0, 3.0]];
        let charges = vec![0.5];
        let types = vec![1u32];
        let mol_ids = vec![1u32];
        let s = LammpsDataWriter::write_atoms_full(&pos, &charges, &types, &mol_ids);
        assert!(s.contains("Atoms"));
        assert!(s.contains("full"));
        assert!(s.contains("0.5"));
    }
    #[test]
    fn test_write_atoms_molecular() {
        let pos = vec![[1.0, 2.0, 3.0]];
        let types = vec![1u32];
        let mol_ids = vec![2u32];
        let s = LammpsDataWriter::write_atoms_molecular(&pos, &types, &mol_ids);
        assert!(s.contains("molecular"));
        assert!(s.contains("2"));
    }
    #[test]
    fn test_lj_cut() {
        let s = LammpsPairStyle::lj_cut(12.0);
        assert!(s.contains("pair_style lj/cut 12"));
    }
    #[test]
    fn test_pair_coeff_lj() {
        let s = LammpsPairStyle::pair_coeff_lj(1, 2, 0.1, 3.4);
        assert!(s.contains("pair_coeff 1 2"));
        assert!(s.contains("0.1"));
        assert!(s.contains("3.4"));
    }
    #[test]
    fn test_lj_cut_coul_long() {
        let s = LammpsPairStyle::lj_cut_coul_long(10.0, 12.0);
        assert!(s.contains("lj/cut/coul/long"));
    }
    #[test]
    fn test_pair_hybrid() {
        let s = LammpsPairStyle::hybrid(&["lj/cut 10.0", "coul/long 12.0"]);
        assert!(s.contains("hybrid"));
    }
    #[test]
    fn test_fix_nve() {
        let s = LammpsFix::nve("1", "all");
        assert!(s.contains("fix 1 all nve"));
    }
    #[test]
    fn test_fix_nvt() {
        let s = LammpsFix::nvt("1", "all", 300.0, 300.0, 100.0);
        assert!(s.contains("nvt"));
        assert!(s.contains("300"));
    }
    #[test]
    fn test_fix_langevin() {
        let s = LammpsFix::langevin("therm", "all", 300.0, 300.0, 100.0, 12345);
        assert!(s.contains("langevin"));
        assert!(s.contains("12345"));
    }
    #[test]
    fn test_fix_setforce() {
        let s = LammpsFix::setforce("f1", "group1", 0.0, 0.0, -1.0);
        assert!(s.contains("setforce"));
        assert!(s.contains("-1"));
    }
    #[test]
    fn test_fix_rigid() {
        let s = LammpsFix::rigid("r1", "all", "single");
        assert!(s.contains("rigid single"));
    }
    #[test]
    fn test_compute_temp() {
        let s = LammpsCompute::temp("mytemp", "all");
        assert!(s.contains("compute mytemp all temp"));
    }
    #[test]
    fn test_compute_pe() {
        let s = LammpsCompute::pe("mype", "all");
        assert!(s.contains("pe"));
    }
    #[test]
    fn test_compute_ke() {
        let s = LammpsCompute::ke("myke", "all");
        assert!(s.contains("ke"));
    }
    #[test]
    fn test_compute_pressure() {
        let s = LammpsCompute::pressure("mypres", "all", "mytemp");
        assert!(s.contains("pressure"));
        assert!(s.contains("mytemp"));
    }
    #[test]
    fn test_compute_rdf() {
        let s = LammpsCompute::rdf("myrdf", "all", 100, 5.0);
        assert!(s.contains("rdf"));
        assert!(s.contains("100"));
    }
    #[test]
    fn test_nve_script() {
        let s = LammpsInputScript::nve_script(5000, 0.001);
        assert!(s.contains("nve"));
        assert!(s.contains("5000"));
        assert!(s.contains("0.001"));
    }
    #[test]
    fn test_thermalise_script() {
        let s = LammpsInputScript::thermalise_script(300.0, 100.0, 10000, 42);
        assert!(s.contains("langevin"));
        assert!(s.contains("300"));
        assert!(s.contains("42"));
    }
    #[test]
    fn test_dump_command() {
        let s = LammpsInputScript::dump_command(
            "d1",
            "all",
            "custom",
            100,
            "dump.lammpstrj",
            "id type x y z",
        );
        assert!(s.contains("dump d1 all custom 100"));
    }
    #[test]
    fn test_thermo_output() {
        let s = LammpsInputScript::thermo_output(100, &["step", "temp", "pe"]);
        assert!(s.contains("thermo 100"));
        assert!(s.contains("thermo_style custom step temp pe"));
    }
    #[test]
    fn test_units_command() {
        let s = LammpsInputScript::units("real");
        assert!(s.contains("units real"));
    }
    #[test]
    fn test_boundary_command() {
        let s = LammpsInputScript::boundary("p p p");
        assert!(s.contains("boundary p p p"));
    }
    #[test]
    fn test_read_data_command() {
        let s = LammpsInputScript::read_data("system.data");
        assert!(s.contains("read_data system.data"));
    }
    #[test]
    fn test_lj_pair_potential_minimum() {
        let pair = LammpsLJPair {
            type_i: 1,
            type_j: 1,
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: None,
        };
        let r_min = pair.r_min();
        let u = pair.potential_energy(r_min);
        assert!((u - (-1.0)).abs() < 1e-10, "U(r_min) should be -ε, got {u}");
    }
    #[test]
    fn test_lj_pair_r_min() {
        let pair = LammpsLJPair {
            type_i: 1,
            type_j: 1,
            epsilon: 0.5,
            sigma: 3.4,
            cutoff: None,
        };
        let r_min = pair.r_min();
        assert!((r_min - 2.0_f64.powf(1.0 / 6.0) * 3.4).abs() < 1e-10);
    }
    #[test]
    fn test_lj_pair_coeff_line() {
        let pair = LammpsLJPair {
            type_i: 1,
            type_j: 2,
            epsilon: 0.1,
            sigma: 3.5,
            cutoff: Some(10.0),
        };
        let s = pair.to_pair_coeff_line();
        assert!(s.contains("pair_coeff 1 2"));
        assert!(s.contains("0.1"));
        assert!(s.contains("3.5"));
        assert!(s.contains("10"));
    }
    #[test]
    fn test_lj_table_lorentz_berthelot() {
        let mut table = LammpsLJTable::new(12.0);
        table.add_pair(LammpsLJPair {
            type_i: 1,
            type_j: 1,
            epsilon: 1.0,
            sigma: 3.0,
            cutoff: None,
        });
        table.add_pair(LammpsLJPair {
            type_i: 2,
            type_j: 2,
            epsilon: 4.0,
            sigma: 5.0,
            cutoff: None,
        });
        table.apply_lorentz_berthelot_mixing();
        let cross = table.get(1, 2).expect("cross-term should be generated");
        assert!(
            (cross.epsilon - 2.0).abs() < 1e-10,
            "ε_12={}",
            cross.epsilon
        );
        assert!((cross.sigma - 4.0).abs() < 1e-10, "σ_12={}", cross.sigma);
    }
    #[test]
    fn test_lj_table_to_input() {
        let mut table = LammpsLJTable::new(10.0);
        table.add_pair(LammpsLJPair {
            type_i: 1,
            type_j: 1,
            epsilon: 0.5,
            sigma: 3.4,
            cutoff: None,
        });
        let s = table.to_lammps_input();
        assert!(s.contains("pair_style lj/cut 10"));
        assert!(s.contains("pair_coeff 1 1"));
    }
    #[test]
    fn test_harmonic_bond_energy() {
        let bond = LammpsHarmonicBond {
            bond_type: 1,
            k: 100.0,
            r0: 1.5,
        };
        assert!((bond.potential_energy(1.5)).abs() < 1e-10);
        assert!((bond.potential_energy(1.6) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_harmonic_bond_coeff_line() {
        let bond = LammpsHarmonicBond {
            bond_type: 2,
            k: 200.0,
            r0: 1.2,
        };
        let s = bond.to_bond_coeff_line();
        assert!(s.contains("bond_coeff 2 200 1.2"));
    }
    #[test]
    fn test_harmonic_angle_energy() {
        let angle = LammpsHarmonicAngle {
            angle_type: 1,
            k: 50.0,
            theta0_deg: 120.0,
        };
        assert!(angle.potential_energy_deg(120.0).abs() < 1e-10);
    }
    #[test]
    fn test_harmonic_angle_coeff_line() {
        let angle = LammpsHarmonicAngle {
            angle_type: 1,
            k: 50.0,
            theta0_deg: 109.47,
        };
        let s = angle.to_angle_coeff_line();
        assert!(s.contains("angle_coeff 1 50 109.47"));
    }
    #[test]
    fn test_log_parser_thermo() {
        let log = "\
LAMMPS (2022)
Step Temp E_pair E_mol TotEng Press
0 300.0 -1000.0 0.0 -800.0 1.0
100 310.5 -999.0 0.0 -798.5 1.1
Loop time
";
        let records = LammpsLogParser::parse_thermo(log);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].step, 0);
        assert_eq!(records[1].step, 100);
    }
    #[test]
    fn test_log_parser_thermo_standard_keywords() {
        let log = "\
Step Temp PE KE Etotal Press Vol
0 300.0 -500.0 100.0 -400.0 1.0 1000.0
200 305.0 -498.0 102.0 -396.0 1.05 1001.0
";
        let records = LammpsLogParser::parse_thermo(log);
        assert_eq!(records.len(), 2);
        assert!((records[0].temp.unwrap() - 300.0).abs() < 1e-6);
        assert!((records[1].vol.unwrap() - 1001.0).abs() < 1e-6);
    }
    #[test]
    fn test_log_parser_count_steps() {
        let log = "\
run 10000
run 5000
minimize 1e-4 1e-8 100 1000
run 20000
";
        let steps = LammpsLogParser::count_total_steps(log);
        assert_eq!(steps, 35000);
    }
    #[test]
    fn test_log_parser_wall_time() {
        let log = "\
LAMMPS run complete
Total wall time: 0:02:30
";
        let wt = LammpsLogParser::parse_wall_time(log);
        assert!(wt.is_some());
        assert!(wt.unwrap().contains("0:02:30"));
    }
    #[test]
    fn test_log_parser_no_wall_time() {
        let wt = LammpsLogParser::parse_wall_time("No time here\n");
        assert!(wt.is_none());
    }
    #[test]
    fn test_lammps_box_volume() {
        let b = LammpsBox::orthogonal([0.0; 3], [10.0, 5.0, 2.0]);
        assert!((b.volume() - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_lammps_box_number_density() {
        let b = LammpsBox::orthogonal([0.0; 3], [10.0, 10.0, 10.0]);
        let rho = b.number_density(1000);
        assert!((rho - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_lammps_box_contains() {
        let b = LammpsBox::orthogonal([0.0; 3], [10.0; 3]);
        assert!(b.contains([5.0, 5.0, 5.0]));
        assert!(!b.contains([11.0, 5.0, 5.0]));
        assert!(!b.contains([-1.0, 5.0, 5.0]));
    }
    #[test]
    fn test_lammps_box_wrap_pbc() {
        let b = LammpsBox::orthogonal([0.0; 3], [10.0; 3]);
        let wrapped = b.wrap_pbc([11.5, -0.5, 5.0]);
        assert!((wrapped[0] - 1.5).abs() < 1e-10, "x={}", wrapped[0]);
        assert!((wrapped[1] - 9.5).abs() < 1e-10, "y={}", wrapped[1]);
        assert!((wrapped[2] - 5.0).abs() < 1e-10, "z={}", wrapped[2]);
    }
    #[test]
    fn test_lammps_box_region_command() {
        let b = LammpsBox::orthogonal([-5.0; 3], [5.0; 3]);
        let s = b.to_region_command("mybox");
        assert!(s.contains("region mybox block"));
        assert!(s.contains("-5"));
        assert!(s.contains("5"));
    }
    #[test]
    fn test_lammps_box_lengths() {
        let b = LammpsBox::orthogonal([1.0, 2.0, 3.0], [4.0, 7.0, 13.0]);
        let l = b.lengths();
        assert!((l[0] - 3.0).abs() < 1e-10);
        assert!((l[1] - 5.0).abs() < 1e-10);
        assert!((l[2] - 10.0).abs() < 1e-10);
    }
}
pub(super) fn read_i32_le(data: &[u8], pos: &mut usize) -> Result<i32> {
    if *pos + 4 > data.len() {
        return Err(Error::Parse("EOF reading i32".into()));
    }
    let v = i32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
    *pos += 4;
    Ok(v)
}
pub(super) fn read_i64_le(data: &[u8], pos: &mut usize) -> Result<i64> {
    if *pos + 8 > data.len() {
        return Err(Error::Parse("EOF reading i64".into()));
    }
    let b = &data[*pos..*pos + 8];
    let v = i64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
    *pos += 8;
    Ok(v)
}
pub(super) fn read_f64_le(data: &[u8], pos: &mut usize) -> Result<f64> {
    if *pos + 8 > data.len() {
        return Err(Error::Parse("EOF reading f64".into()));
    }
    let b = &data[*pos..*pos + 8];
    let v = f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
    *pos += 8;
    Ok(v)
}
/// Parse a line in `atom_style full` format used by OVITO:
/// `atom_id mol_id atom_type charge x y z`
///
/// Returns `(atom_id, mol_id, atom_type, charge, x, y, z)` on success.
pub fn parse_atom_style_full(line: &str) -> Option<(usize, usize, usize, f64, f64, f64, f64)> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 7 {
        return None;
    }
    let atom_id = fields[0].parse::<usize>().ok()?;
    let mol_id = fields[1].parse::<usize>().ok()?;
    let atom_type = fields[2].parse::<usize>().ok()?;
    let charge = fields[3].parse::<f64>().ok()?;
    let x = fields[4].parse::<f64>().ok()?;
    let y = fields[5].parse::<f64>().ok()?;
    let z = fields[6].parse::<f64>().ok()?;
    Some((atom_id, mol_id, atom_type, charge, x, y, z))
}
/// Parse a line in `atom_style charge` format:
/// `atom_id atom_type charge x y z`
pub fn parse_atom_style_charge(line: &str) -> Option<(usize, usize, f64, f64, f64, f64)> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 6 {
        return None;
    }
    let atom_id = fields[0].parse::<usize>().ok()?;
    let atom_type = fields[1].parse::<usize>().ok()?;
    let charge = fields[2].parse::<f64>().ok()?;
    let x = fields[3].parse::<f64>().ok()?;
    let y = fields[4].parse::<f64>().ok()?;
    let z = fields[5].parse::<f64>().ok()?;
    Some((atom_id, atom_type, charge, x, y, z))
}
#[cfg(test)]
mod tests_lammps_additions {
    use super::*;
    use crate::LammpsAtom;

    use crate::lammps::DumpColumnSpec;
    use crate::lammps::LammpsBinaryDumpReader;
    use crate::lammps::LammpsBinaryDumpWriter;
    use crate::lammps::LammpsDataSectionBuilder;
    use crate::lammps::LammpsRestartReader;
    use crate::lammps::LammpsRestartWriter;
    use crate::lammps::LammpsRunBlock;

    use oxiphysics_core::math::Vec3;
    fn make_atom(id: usize, x: f64, y: f64, z: f64) -> LammpsAtom {
        LammpsAtom {
            id,
            type_id: 1,
            position: Vec3::new(x, y, z),
            velocity: Vec3::new(0.0, 0.0, 0.0),
        }
    }
    #[test]
    fn test_binary_dump_frame_size() {
        assert_eq!(
            LammpsBinaryDumpWriter::frame_byte_size(2),
            8 + 8 + 48 + 2 * 56
        );
    }
    #[test]
    fn test_binary_dump_round_trip() {
        let atoms = vec![make_atom(1, 1.0, 2.0, 3.0), make_atom(2, 4.0, 5.0, 6.0)];
        let box_bounds = [[0.0, 10.0]; 3];
        let mut buf: Vec<u8> = Vec::new();
        LammpsBinaryDumpWriter::write_frame(&mut buf, 42, box_bounds, &atoms).unwrap();
        let (ts, bb, read_atoms, _) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
        assert_eq!(ts, 42);
        assert_eq!(read_atoms.len(), 2);
        assert!((read_atoms[0].position.x - 1.0).abs() < 1e-12);
        assert!((read_atoms[1].position.z - 6.0).abs() < 1e-12);
        assert!((bb[0][0] - 0.0).abs() < 1e-12);
        assert!((bb[2][1] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_binary_dump_empty_atoms() {
        let mut buf: Vec<u8> = Vec::new();
        LammpsBinaryDumpWriter::write_frame(&mut buf, 0, [[0.0, 1.0]; 3], &[]).unwrap();
        let (ts, _bb, atoms, _) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
        assert_eq!(ts, 0);
        assert!(atoms.is_empty());
    }
    #[test]
    fn test_binary_dump_truncated_returns_error() {
        let buf = vec![0u8; 4];
        assert!(LammpsBinaryDumpReader::read_frame(&buf, 0).is_err());
    }
    #[test]
    fn test_binary_dump_ids_preserved() {
        let atoms = vec![make_atom(99, 0.0, 0.0, 0.0), make_atom(200, 1.0, 1.0, 1.0)];
        let mut buf: Vec<u8> = Vec::new();
        LammpsBinaryDumpWriter::write_frame(&mut buf, 1, [[0.0, 5.0]; 3], &atoms).unwrap();
        let (_, _, read_atoms, _) = LammpsBinaryDumpReader::read_frame(&buf, 0).unwrap();
        assert_eq!(read_atoms[0].id, 99);
        assert_eq!(read_atoms[1].id, 200);
    }
    #[test]
    fn test_restart_round_trip() {
        let atoms = vec![make_atom(1, 1.5, 2.5, 3.5)];
        let masses = vec![12.0, 16.0];
        let mut buf: Vec<u8> = Vec::new();
        LammpsRestartWriter::write(&mut buf, 1000, &masses, &atoms).unwrap();
        let (ts, m, a) = LammpsRestartReader::read(&buf).unwrap();
        assert_eq!(ts, 1000);
        assert_eq!(m.len(), 2);
        assert!((m[0] - 12.0).abs() < 1e-12);
        assert_eq!(a.len(), 1);
        assert!((a[0].position.y - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_restart_bad_magic() {
        let buf = b"XRST\x01\x00\x00\x00".to_vec();
        assert!(LammpsRestartReader::read(&buf).is_err());
    }
    #[test]
    fn test_restart_empty_atoms() {
        let mut buf: Vec<u8> = Vec::new();
        LammpsRestartWriter::write(&mut buf, 5, &[1.0], &[]).unwrap();
        let (ts, masses, atoms) = LammpsRestartReader::read(&buf).unwrap();
        assert_eq!(ts, 5);
        assert_eq!(masses.len(), 1);
        assert!(atoms.is_empty());
    }
    #[test]
    fn test_dump_column_spec_standard() {
        let spec = DumpColumnSpec::standard();
        assert_eq!(spec.len(), 8);
        assert!(spec.contains("id"));
        assert!(spec.contains("vz"));
        assert!(spec.header_line().starts_with("ITEM: ATOMS"));
    }
    #[test]
    fn test_dump_column_spec_positions_only() {
        let spec = DumpColumnSpec::positions_only();
        assert_eq!(spec.len(), 5);
        assert!(!spec.contains("vx"));
        assert!(spec.contains("z"));
    }
    #[test]
    fn test_dump_column_spec_index_of() {
        let spec = DumpColumnSpec::standard();
        assert_eq!(spec.index_of("x"), Some(2));
        assert_eq!(spec.index_of("vz"), Some(7));
        assert_eq!(spec.index_of("missing"), None);
    }
    #[test]
    fn test_dump_column_spec_header_line() {
        let spec = DumpColumnSpec::standard();
        let h = spec.header_line();
        assert_eq!(h, "ITEM: ATOMS id type x y z vx vy vz");
    }
    #[test]
    fn test_data_section_masses() {
        let mut b = LammpsDataSectionBuilder::new();
        b.add_mass(1, 12.011);
        b.add_mass(2, 15.999);
        let s = b.masses_section();
        assert!(s.starts_with("Masses"));
        assert!(s.contains("12.011000"));
        assert!(s.contains("15.999000"));
    }
    #[test]
    fn test_data_section_pair_coeffs() {
        let mut b = LammpsDataSectionBuilder::new();
        b.add_pair_coeff(1, 1, 0.086, 3.4);
        let s = b.pair_coeffs_section();
        assert!(s.contains("Pair Coeffs"));
        assert!(s.contains("0.086000"));
        assert!(s.contains("3.400000"));
    }
    #[test]
    fn test_data_section_bond_coeffs() {
        let mut b = LammpsDataSectionBuilder::new();
        b.add_bond_coeff(1, 350.0, 1.5);
        let s = b.bond_coeffs_section();
        assert!(s.contains("Bond Coeffs"));
        assert!(s.contains("350.000000"));
        assert!(s.contains("1.500000"));
    }
    #[test]
    fn test_data_section_build_all() {
        let mut b = LammpsDataSectionBuilder::new();
        b.add_mass(1, 12.0);
        b.add_pair_coeff(1, 1, 0.1, 3.0);
        b.add_bond_coeff(1, 100.0, 1.4);
        b.add_angle_coeff(1, 50.0, 109.5);
        let s = b.build();
        assert!(s.contains("Masses"));
        assert!(s.contains("Pair Coeffs"));
        assert!(s.contains("Bond Coeffs"));
        assert!(s.contains("Angle Coeffs"));
    }
    #[test]
    fn test_parse_atom_style_full() {
        let line = "1 1 2 -0.5 1.2 3.4 5.6";
        let r = parse_atom_style_full(line).unwrap();
        assert_eq!(r.0, 1);
        assert_eq!(r.1, 1);
        assert_eq!(r.2, 2);
        assert!((r.3 - (-0.5)).abs() < 1e-12);
        assert!((r.4 - 1.2).abs() < 1e-12);
    }
    #[test]
    fn test_parse_atom_style_full_too_short() {
        assert!(parse_atom_style_full("1 2 3").is_none());
    }
    #[test]
    fn test_parse_atom_style_charge() {
        let line = "5 3 1.0 0.0 2.0 4.0";
        let r = parse_atom_style_charge(line).unwrap();
        assert_eq!(r.0, 5);
        assert_eq!(r.1, 3);
        assert!((r.2 - 1.0).abs() < 1e-12);
        assert!((r.5 - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_parse_atom_style_charge_bad() {
        assert!(parse_atom_style_charge("abc def ghi jkl mno").is_none());
    }
    #[test]
    fn test_run_block_nve_script() {
        let block = LammpsRunBlock::nve(1.0, 10000);
        let s = block.to_script();
        assert!(s.contains("timestep 1.000000"), "should set timestep");
        assert!(s.contains("fix 1 all nve"), "should set nve fix");
        assert!(s.contains("run 10000"), "should contain run");
    }
    #[test]
    fn test_run_block_nvt_script() {
        let block = LammpsRunBlock::nvt(2.0, 5000, 300.0);
        let s = block.to_script();
        assert!(s.contains("nvt"), "should contain nvt");
        assert!(s.contains("300"), "should contain temperature");
    }
    #[test]
    fn test_run_block_npt_script() {
        let block = LammpsRunBlock::npt(2.0, 5000, 300.0, 1.0);
        let s = block.to_script();
        assert!(s.contains("npt"), "should contain npt");
        assert!(s.contains("300"), "should contain temperature");
        assert!(
            s.contains("1.0") || s.contains("1."),
            "should contain pressure"
        );
    }
    #[test]
    fn test_run_block_total_time() {
        let block = LammpsRunBlock::nve(0.5, 2000);
        assert!((block.total_time_fs() - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_run_block_custom_freqs() {
        let block = LammpsRunBlock::nvt(1.0, 1000, 298.0)
            .thermo_every(50)
            .dump_every(200);
        assert_eq!(block.thermo_freq, 50);
        assert_eq!(block.dump_freq, 200);
        let s = block.to_script();
        assert!(s.contains("thermo 50"), "should have thermo freq");
        assert!(s.contains("200"), "dump line should mention dump freq");
    }
}
