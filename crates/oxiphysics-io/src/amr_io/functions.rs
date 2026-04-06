//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use super::types::{AmrHierarchy, AmrWriter};

pub(super) fn parse_attr<'a>(s: &'a str, attr: &str) -> Option<&'a str> {
    let needle = format!(r#"{attr}=""#);
    let start = s.find(needle.as_str())? + needle.len();
    let rest = &s[start..];
    let end = rest.find('"')?;
    Some(&rest[..end])
}
pub(super) fn parse_inline_value(s: &str) -> Option<&str> {
    let start = s.find('>')? + 1;
    let end = s[start..].find('<').map(|p| start + p)?;
    if start < end {
        Some(&s[start..end])
    } else {
        None
    }
}
/// Write a PVD collection file pointing to per-timestep VTK files.
///
/// `path` — path to the .pvd file to create.
/// `hierarchy` — hierarchy used for metadata (time/cycle embedded in filename).
pub fn write_amr_pvd(path: &str, hierarchy: &AmrHierarchy) -> io::Result<()> {
    let dir = Path::new(path).parent().unwrap_or_else(|| Path::new("."));
    let stem = Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("amr");
    let vtk_name = format!("{stem}_t{:.6}_c{}.vthb", hierarchy.time, hierarchy.cycle);
    let vtk_path = dir.join(&vtk_name);
    AmrWriter::write_vtk(vtk_path.to_str().unwrap_or(&vtk_name), hierarchy)?;
    let mut f = File::create(path)?;
    writeln!(f, r#"<?xml version="1.0"?>"#)?;
    writeln!(
        f,
        r#"<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">"#
    )?;
    writeln!(f, r#"  <Collection>"#)?;
    writeln!(
        f,
        r#"    <DataSet timestep="{:.6}" group="" part="0" file="{}"/>"#,
        hierarchy.time, vtk_name
    )?;
    writeln!(f, r#"  </Collection>"#)?;
    writeln!(f, r#"</VTKFile>"#)?;
    Ok(())
}
pub(super) const CHECKPOINT_MAGIC: u32 = 0x414D_5243;
pub(super) const CHECKPOINT_VERSION: u32 = 1;
/// Morton code computed from 3-D integer coordinates.
///
/// Interleaves the bits of `(i, j, k)` using up to 21 bits per axis to
/// produce a 63-bit Z-curve key.
pub fn morton_encode(i: u32, j: u32, k: u32) -> u64 {
    fn expand(v: u32) -> u64 {
        let mut x = v as u64 & 0x001f_ffff;
        x = (x | (x << 32)) & 0x001f_0000_0000_ffff;
        x = (x | (x << 16)) & 0x001f_0000_ff00_00ff;
        x = (x | (x << 8)) & 0x100f_00f0_0f00_f00f;
        x = (x | (x << 4)) & 0x10c3_0c30_c30c_30c3;
        x = (x | (x << 2)) & 0x1249_2492_4924_9249;
        x
    }
    expand(i) | (expand(j) << 1) | (expand(k) << 2)
}
pub(super) fn extract_f64_after(s: &str, key: &str) -> Option<f64> {
    let pos = s.find(key)? + key.len();
    let rest = s[pos..].trim_start();
    let end = rest
        .find(|c: char| {
            !c.is_ascii_digit() && c != '.' && c != '-' && c != 'e' && c != 'E' && c != '+'
        })
        .unwrap_or(rest.len());
    rest[..end].parse::<f64>().ok()
}
pub(super) fn extract_u64_after(s: &str, key: &str) -> Option<u64> {
    let pos = s.find(key)? + key.len();
    let rest = s[pos..].trim_start();
    let end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    rest[..end].parse::<u64>().ok()
}
pub(super) fn extract_usize_at(s: &str) -> Option<usize> {
    let rest = s.trim_start();
    let end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    rest[..end].parse::<usize>().ok()
}
pub(super) fn extract_f64_array_6(s: &str) -> Option<[f64; 6]> {
    let start = s.find('[')? + 1;
    let end = s.find(']')?;
    let nums: Vec<f64> = s[start..end]
        .split(',')
        .filter_map(|t| t.trim().parse::<f64>().ok())
        .collect();
    if nums.len() < 6 {
        return None;
    }
    Some([nums[0], nums[1], nums[2], nums[3], nums[4], nums[5]])
}
pub(super) fn extract_i64_array_3(s: &str) -> Option<[i64; 3]> {
    let key = "\"coords\":";
    let pos = s.find(key)? + key.len();
    let start = s[pos..].find('[')? + pos + 1;
    let end = s[start..].find(']')? + start;
    let nums: Vec<i64> = s[start..end]
        .split(',')
        .filter_map(|t| t.trim().parse::<i64>().ok())
        .collect();
    if nums.len() < 3 {
        return None;
    }
    Some([nums[0], nums[1], nums[2]])
}
pub(super) fn extract_f64_array_vec(s: &str) -> Vec<f64> {
    let key = "\"data\":";
    let pos = match s.find(key) {
        Some(p) => p + key.len(),
        None => return Vec::new(),
    };
    let start = match s[pos..].find('[') {
        Some(p) => pos + p + 1,
        None => return Vec::new(),
    };
    let end = match s[start..].find(']') {
        Some(p) => start + p,
        None => return Vec::new(),
    };
    s[start..end]
        .split(',')
        .filter_map(|t| t.trim().parse::<f64>().ok())
        .collect()
}
pub(super) fn extract_string_array(s: &str) -> Vec<String> {
    let start = match s.find('[') {
        Some(p) => p + 1,
        None => return Vec::new(),
    };
    let end = match s.find(']') {
        Some(p) => p,
        None => return Vec::new(),
    };
    let inner = &s[start..end];
    let mut result = Vec::new();
    let mut from = 0;
    while let Some(q1) = inner[from..].find('"') {
        let abs_q1 = from + q1 + 1;
        if let Some(q2) = inner[abs_q1..].find('"') {
            result.push(inner[abs_q1..abs_q1 + q2].to_string());
            from = abs_q1 + q2 + 1;
        } else {
            break;
        }
    }
    result
}
pub(super) fn find_matching_bracket(s: &str) -> Option<usize> {
    let mut depth = 0i32;
    let mut found_open = false;
    for (i, c) in s.char_indices() {
        match c {
            '[' => {
                depth += 1;
                found_open = true;
            }
            ']' => {
                depth -= 1;
                if found_open && depth == 0 {
                    return Some(i);
                }
            }
            _ => {}
        }
    }
    None
}
pub(super) fn find_matching_brace(s: &str) -> Option<usize> {
    let mut depth = 0i32;
    let mut found_open = false;
    for (i, c) in s.char_indices() {
        match c {
            '{' => {
                depth += 1;
                found_open = true;
            }
            '}' => {
                depth -= 1;
                if found_open && depth == 0 {
                    return Some(i);
                }
            }
            _ => {}
        }
    }
    None
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::amr_io::AmrBlock;
    use crate::amr_io::AmrCell;
    use crate::amr_io::AmrCheckpoint;
    use crate::amr_io::AmrCoarsening;
    use crate::amr_io::AmrFieldInterp;
    use crate::amr_io::AmrGrid;
    use crate::amr_io::AmrLevel;
    use crate::amr_io::AmrMeshExport;
    use crate::amr_io::AmrReader;
    use crate::amr_io::AmrRefinementCriteria;
    use crate::amr_io::AmrSerializer;
    use crate::amr_io::AmrStats;
    use crate::amr_io::AmrtNode;
    use crate::amr_io::AmrtTree;
    use crate::amr_io::RefinementMethod;
    use crate::amr_io::types::*;
    use std::io::Read;
    fn make_two_level_hierarchy() -> AmrHierarchy {
        let mut h = AmrHierarchy::new();
        h.time = 1.5;
        h.cycle = 42;
        let lm0 = AmrLevel::new(0, 1.0, [0.0, 0.0, 0.0, 4.0, 4.0, 4.0]);
        let mut g0 = AmrGrid::new(lm0, vec!["pressure".into(), "density".into()]);
        g0.add_cell(AmrCell::new(0, [0, 0, 0], vec![1.0, 2.0]));
        g0.add_cell(AmrCell::new(0, [1, 0, 0], vec![3.0, 4.0]));
        g0.add_cell(AmrCell::new(0, [0, 1, 0], vec![5.0, 6.0]));
        h.add_level(g0);
        let lm1 = AmrLevel::new(1, 0.5, [0.0, 0.0, 0.0, 4.0, 4.0, 4.0]);
        let mut g1 = AmrGrid::new(lm1, vec!["pressure".into(), "density".into()]);
        g1.add_cell(AmrCell::new(1, [0, 0, 0], vec![1.1, 2.1]));
        g1.add_cell(AmrCell::new(1, [1, 0, 0], vec![1.2, 2.2]));
        g1.add_cell(AmrCell::new(1, [0, 1, 0], vec![1.3, 2.3]));
        g1.add_cell(AmrCell::new(1, [1, 1, 0], vec![1.4, 2.4]));
        h.add_level(g1);
        h
    }
    #[test]
    fn test_amr_level_new() {
        let lm = AmrLevel::new(2, 0.25, [0.0; 6]);
        assert_eq!(lm.level, 2);
        assert!((lm.cell_size - 0.25).abs() < 1e-12);
    }
    #[test]
    fn test_refinement_ratio() {
        let lm0 = AmrLevel::new(0, 1.0, [0.0; 6]);
        let lm2 = AmrLevel::new(2, 0.25, [0.0; 6]);
        assert_eq!(lm0.refinement_ratio(), 1);
        assert_eq!(lm2.refinement_ratio(), 4);
    }
    #[test]
    fn test_amr_level_domain() {
        let bounds = [0.0, 0.0, 0.0, 4.0, 4.0, 1.0];
        let lm = AmrLevel::new(0, 1.0, bounds);
        assert!((lm.domain_bounds[3] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_amr_cell_new() {
        let c = AmrCell::new(0, [1, 2, 3], vec![1.0, 2.0]);
        assert_eq!(c.coords, [1, 2, 3]);
        assert_eq!(c.data.len(), 2);
    }
    #[test]
    fn test_physical_center() {
        let lm = AmrLevel::new(0, 1.0, [0.0, 0.0, 0.0, 4.0, 4.0, 4.0]);
        let c = AmrCell::new(0, [0, 0, 0], vec![]);
        let center = c.physical_center(&lm);
        assert!((center[0] - 0.5).abs() < 1e-12);
        assert!((center[1] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_physical_center_offset() {
        let lm = AmrLevel::new(0, 0.5, [1.0, 2.0, 3.0, 5.0, 6.0, 7.0]);
        let c = AmrCell::new(0, [2, 1, 0], vec![]);
        let center = c.physical_center(&lm);
        assert!((center[0] - (1.0 + 2.5 * 0.5)).abs() < 1e-10);
    }
    #[test]
    fn test_amr_grid_cell_count() {
        let h = make_two_level_hierarchy();
        assert_eq!(h.levels[0].cell_count(), 3);
        assert_eq!(h.levels[1].cell_count(), 4);
    }
    #[test]
    fn test_amr_grid_bounding_box() {
        let h = make_two_level_hierarchy();
        let bb = h.levels[0].integer_bounding_box().unwrap();
        assert_eq!(bb[0], 0);
        assert_eq!(bb[3], 1);
        assert_eq!(bb[4], 1);
    }
    #[test]
    fn test_amr_grid_physical_bounding_box() {
        let h = make_two_level_hierarchy();
        let pbb = h.levels[0].physical_bounding_box().unwrap();
        assert!((pbb[0] - 0.0).abs() < 1e-12);
        assert!((pbb[3] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_amr_grid_empty_bounding_box() {
        let lm = AmrLevel::new(0, 1.0, [0.0; 6]);
        let grid = AmrGrid::new(lm, vec![]);
        assert!(grid.integer_bounding_box().is_none());
    }
    #[test]
    fn test_amr_grid_field_names() {
        let h = make_two_level_hierarchy();
        assert_eq!(h.levels[0].field_names[0], "pressure");
        assert_eq!(h.levels[0].field_names[1], "density");
    }
    #[test]
    fn test_hierarchy_num_levels() {
        let h = make_two_level_hierarchy();
        assert_eq!(h.num_levels(), 2);
    }
    #[test]
    fn test_hierarchy_total_cells() {
        let h = make_two_level_hierarchy();
        assert_eq!(h.total_cells(), 7);
    }
    #[test]
    fn test_hierarchy_time_cycle() {
        let h = make_two_level_hierarchy();
        assert!((h.time - 1.5).abs() < 1e-12);
        assert_eq!(h.cycle, 42);
    }
    #[test]
    fn test_coarsen_field_basic() {
        let h = make_two_level_hierarchy();
        let coarsened = h.coarsen_field(1, 0);
        assert_eq!(coarsened.len(), 3);
        let expected = (1.1 + 1.2 + 1.3 + 1.4) / 4.0;
        assert!((coarsened[0] - expected).abs() < 1e-10);
    }
    #[test]
    fn test_refine_field_basic() {
        let h = make_two_level_hierarchy();
        let refined = h.refine_field(0, 0);
        assert_eq!(refined.len(), 4);
        assert!((refined[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_coarsen_out_of_bounds() {
        let h = make_two_level_hierarchy();
        let r = h.coarsen_field(0, 0);
        assert!(r.is_empty());
    }
    #[test]
    fn test_refine_out_of_bounds() {
        let h = make_two_level_hierarchy();
        let r = h.refine_field(1, 0);
        assert!(r.is_empty());
    }
    #[test]
    fn test_amr_stats_cells_per_level() {
        let h = make_two_level_hierarchy();
        let stats = AmrStats::compute(&h);
        assert_eq!(stats.cells_per_level, vec![3, 4]);
    }
    #[test]
    fn test_amr_stats_total() {
        let h = make_two_level_hierarchy();
        let stats = AmrStats::compute(&h);
        assert_eq!(stats.total_cells, 7);
    }
    #[test]
    fn test_amr_stats_memory() {
        let h = make_two_level_hierarchy();
        let stats = AmrStats::compute(&h);
        assert_eq!(stats.memory_bytes, 7 * 2 * 8);
    }
    #[test]
    fn test_interpolate_at_fine() {
        let h = make_two_level_hierarchy();
        let val = AmrFieldInterp::interpolate_at(&h, [0.25, 0.25, 0.0], 0);
        assert!(val.is_some());
        assert!((val.unwrap() - 1.1).abs() < 1e-10);
    }
    #[test]
    fn test_interpolate_at_miss() {
        let h = make_two_level_hierarchy();
        let val = AmrFieldInterp::interpolate_at(&h, [100.0, 100.0, 0.0], 0);
        assert!(val.is_none());
    }
    #[test]
    fn test_bilinear_xy() {
        let h = make_two_level_hierarchy();
        let grid = &h.levels[1];
        let val = AmrFieldInterp::bilinear_xy(grid, 0.25, 0.25, 0);
        assert!(val.is_some());
    }
    #[test]
    fn test_bilinear_xy_no_cell() {
        let lm = AmrLevel::new(0, 1.0, [0.0; 6]);
        let grid = AmrGrid::new(lm, vec!["p".into()]);
        let val = AmrFieldInterp::bilinear_xy(&grid, 0.5, 0.5, 0);
        assert!(val.is_none());
    }
    #[test]
    fn test_vtk_write_and_read_roundtrip() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_roundtrip.vthb";
        AmrWriter::write_vtk(path, &h).expect("write failed");
        let h2 = AmrReader::read_vtk(path).expect("read failed");
        assert_eq!(h2.num_levels(), h.num_levels());
        assert_eq!(h2.total_cells(), h.total_cells());
    }
    #[test]
    fn test_vtk_write_field_values() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_fields.vthb";
        AmrWriter::write_vtk(path, &h).expect("write failed");
        let h2 = AmrReader::read_vtk(path).expect("read failed");
        let cell = &h2.levels[0].cells[0];
        assert!((cell.data[0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_vtk_write_creates_file() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_exists.vthb";
        AmrWriter::write_vtk(path, &h).expect("write failed");
        assert!(std::path::Path::new(path).exists());
    }
    #[test]
    fn test_write_amr_pvd() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr.pvd";
        write_amr_pvd(path, &h).expect("pvd write failed");
        assert!(std::path::Path::new(path).exists());
    }
    #[test]
    fn test_pvd_content_has_timestep() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_ts.pvd";
        write_amr_pvd(path, &h).expect("pvd write failed");
        let mut raw = String::new();
        File::open(path).unwrap().read_to_string(&mut raw).unwrap();
        assert!(raw.contains("timestep="));
        assert!(raw.contains("1.500000"));
    }
    #[test]
    fn test_checkpoint_write_read_roundtrip() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_checkpoint.amrc";
        AmrCheckpoint::write(path, &h).expect("checkpoint write failed");
        let h2 = AmrCheckpoint::read(path).expect("checkpoint read failed");
        assert_eq!(h2.num_levels(), h.num_levels());
        assert_eq!(h2.total_cells(), h.total_cells());
        assert!((h2.time - h.time).abs() < 1e-12);
        assert_eq!(h2.cycle, h.cycle);
    }
    #[test]
    fn test_checkpoint_field_names_preserved() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_ckpt_fields.amrc";
        AmrCheckpoint::write(path, &h).expect("write");
        let h2 = AmrCheckpoint::read(path).expect("read");
        assert_eq!(h2.levels[0].field_names, vec!["pressure", "density"]);
    }
    #[test]
    fn test_checkpoint_cell_data_preserved() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_ckpt_data.amrc";
        AmrCheckpoint::write(path, &h).expect("write");
        let h2 = AmrCheckpoint::read(path).expect("read");
        let orig = &h.levels[0].cells[0].data;
        let restored = &h2.levels[0].cells[0].data;
        assert_eq!(orig.len(), restored.len());
        for (o, r) in orig.iter().zip(restored.iter()) {
            assert!((o - r).abs() < 1e-12);
        }
    }
    #[test]
    fn test_checkpoint_bad_magic() {
        let path = "/tmp/test_amr_bad_magic.amrc";
        std::fs::write(path, b"BAAD\x00\x00\x00\x00").unwrap();
        let result = AmrCheckpoint::read(path);
        assert!(result.is_err());
    }
    #[test]
    fn test_checkpoint_coords_preserved() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_ckpt_coords.amrc";
        AmrCheckpoint::write(path, &h).expect("write");
        let h2 = AmrCheckpoint::read(path).expect("read");
        assert_eq!(h2.levels[1].cells[1].coords, [1, 0, 0]);
    }
    #[test]
    fn test_empty_hierarchy() {
        let h = AmrHierarchy::new();
        assert_eq!(h.num_levels(), 0);
        assert_eq!(h.total_cells(), 0);
    }
    #[test]
    fn test_hierarchy_default() {
        let h = AmrHierarchy::default();
        assert_eq!(h.num_levels(), 0);
    }
    #[test]
    fn test_stats_empty_hierarchy() {
        let h = AmrHierarchy::new();
        let stats = AmrStats::compute(&h);
        assert_eq!(stats.total_cells, 0);
        assert!(stats.cells_per_level.is_empty());
    }
    #[test]
    fn test_amr_cell_empty_data() {
        let c = AmrCell::new(0, [0; 3], vec![]);
        assert!(c.data.is_empty());
    }
    #[test]
    fn test_three_level_hierarchy() {
        let mut h = AmrHierarchy::new();
        for lv in 0..3 {
            let cs = 1.0 / (1u64 << lv) as f64;
            let lm = AmrLevel::new(lv, cs, [0.0; 6]);
            let mut g = AmrGrid::new(lm, vec!["u".into()]);
            g.add_cell(AmrCell::new(lv, [0, 0, 0], vec![lv as f64]));
            h.add_level(g);
        }
        assert_eq!(h.num_levels(), 3);
        assert_eq!(h.total_cells(), 3);
    }
    #[test]
    fn test_vtk_roundtrip_preserves_cell_count_per_level() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_roundtrip2.vthb";
        AmrWriter::write_vtk(path, &h).expect("write");
        let h2 = AmrReader::read_vtk(path).expect("read");
        for (orig_g, restored_g) in h.levels.iter().zip(h2.levels.iter()) {
            assert_eq!(orig_g.cell_count(), restored_g.cell_count());
        }
    }
    #[test]
    fn test_morton_encode_origin() {
        assert_eq!(morton_encode(0, 0, 0), 0);
    }
    #[test]
    fn test_morton_encode_different_coords() {
        let a = morton_encode(1, 0, 0);
        let b = morton_encode(0, 1, 0);
        assert_ne!(a, b);
    }
    #[test]
    fn test_morton_encode_sorted() {
        let c1 = morton_encode(0, 0, 0);
        let c2 = morton_encode(1, 0, 0);
        assert!(c2 > c1);
    }
    #[test]
    fn test_amrt_tree_init_root() {
        let mut tree = AmrtTree::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 3);
        tree.init_root();
        assert_eq!(tree.nodes.len(), 1);
        assert!(tree.nodes[0].is_leaf);
        assert_eq!(tree.leaf_count(), 1);
    }
    #[test]
    fn test_amrt_tree_refine_root() {
        let mut tree = AmrtTree::new([0.0, 0.0, 0.0, 8.0, 8.0, 8.0], 3);
        tree.init_root();
        let ok = tree.refine(0);
        assert!(ok);
        assert_eq!(tree.leaf_count(), 8);
    }
    #[test]
    fn test_amrt_tree_refine_max_level() {
        let mut tree = AmrtTree::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 0);
        tree.init_root();
        let ok = tree.refine(0);
        assert!(!ok);
    }
    #[test]
    fn test_amrt_tree_coarsen() {
        let mut tree = AmrtTree::new([0.0, 0.0, 0.0, 8.0, 8.0, 8.0], 3);
        tree.init_root();
        tree.refine(0);
        let ok = tree.coarsen(0);
        assert!(ok);
        assert_eq!(tree.leaf_count(), 1);
    }
    #[test]
    fn test_amrt_tree_max_active_level() {
        let mut tree = AmrtTree::new([0.0, 0.0, 0.0, 8.0, 8.0, 8.0], 3);
        tree.init_root();
        tree.refine(0);
        assert_eq!(tree.max_active_level(), 1);
    }
    #[test]
    fn test_amrt_tree_cell_size_at() {
        let tree = AmrtTree::new([0.0, 0.0, 0.0, 8.0, 8.0, 8.0], 3);
        let cs = tree.cell_size_at(0);
        assert!((cs[0] - 8.0_f64).abs() < 1.0e-12_f64);
        let cs1 = tree.cell_size_at(1);
        assert!((cs1[0] - 4.0_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_amrt_tree_node_centre() {
        let mut tree = AmrtTree::new([0.0, 0.0, 0.0, 8.0, 8.0, 8.0], 3);
        tree.init_root();
        let centre = tree.node_centre(0);
        assert!((centre[0] - 4.0_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_amrt_octant_index() {
        assert_eq!(AmrtNode::octant_index(0, 0, 0), 0);
        assert_eq!(AmrtNode::octant_index(1, 0, 0), 1);
        assert_eq!(AmrtNode::octant_index(1, 1, 1), 7);
    }
    #[test]
    fn test_amr_block_new() {
        let b = AmrBlock::new(0, [0, 0, 0], [3, 3, 3], 2);
        assert_eq!(b.dims, [4, 4, 4]);
        assert_eq!(b.num_fields, 2);
        assert_eq!(b.cell_count(), 64);
    }
    #[test]
    fn test_amr_block_set_get() {
        let mut b = AmrBlock::new(0, [0, 0, 0], [1, 1, 1], 1);
        b.set(0, 0, 0, 0, 3.125_f64);
        assert!((b.get(0, 0, 0, 0) - 3.125_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_amr_block_fill_field() {
        let mut b = AmrBlock::new(0, [0, 0, 0], [1, 1, 1], 1);
        b.fill_field(0, 7.0_f64);
        assert!((b.get(1, 1, 1, 0) - 7.0_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_amr_block_field_max_min() {
        let mut b = AmrBlock::new(0, [0, 0, 0], [1, 1, 0], 1);
        b.set(0, 0, 0, 0, 1.0_f64);
        b.set(1, 0, 0, 0, 5.0_f64);
        b.set(0, 1, 0, 0, 3.0_f64);
        b.set(1, 1, 0, 0, 2.0_f64);
        assert!((b.field_max(0) - 5.0_f64).abs() < 1.0e-12_f64);
        assert!((b.field_min(0) - 1.0_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_amr_block_refinement_flag_default_false() {
        let b = AmrBlock::new(1, [0, 0, 0], [2, 2, 2], 1);
        assert!(!b.refinement_flag);
    }
    #[test]
    fn test_amr_serializer_json_roundtrip() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        assert!(json.contains("\"time\""));
        assert!(json.contains("\"cycle\""));
        assert!(json.contains("\"cells\""));
    }
    #[test]
    fn test_amr_serializer_json_contains_time() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        assert!(json.contains("1.5"));
    }
    #[test]
    fn test_amr_serializer_json_contains_field_names() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        assert!(json.contains("pressure"));
        assert!(json.contains("density"));
    }
    #[test]
    fn test_amr_serializer_json_contains_levels() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        assert!(json.contains("\"level\":0"));
        assert!(json.contains("\"level\":1"));
    }
    #[test]
    fn test_amr_serializer_from_json_basic() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        let h2 = AmrSerializer::from_json(&json);
        assert!(h2.is_some());
    }
    #[test]
    fn test_amr_serializer_from_json_time_preserved() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        let h2 = AmrSerializer::from_json(&json).unwrap();
        assert!((h2.time - 1.5_f64).abs() < 1.0e-6_f64);
    }
    #[test]
    fn test_amr_serializer_from_json_cycle_preserved() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        let h2 = AmrSerializer::from_json(&json).unwrap();
        assert_eq!(h2.cycle, 42);
    }
    #[test]
    fn test_amr_serializer_from_json_num_levels() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        let h2 = AmrSerializer::from_json(&json).unwrap();
        assert_eq!(h2.num_levels(), h.num_levels());
    }
    #[test]
    fn test_amr_serializer_from_json_cell_count() {
        let h = make_two_level_hierarchy();
        let json = AmrSerializer::to_json(&h);
        let h2 = AmrSerializer::from_json(&json).unwrap();
        assert_eq!(h2.total_cells(), h.total_cells());
    }
    #[test]
    fn test_amr_mesh_export_creates_file() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_export.vtu";
        AmrMeshExport::write_vtu(path, &h, 0).expect("vtu write failed");
        assert!(std::path::Path::new(path).exists());
    }
    #[test]
    fn test_amr_mesh_export_vtk_header() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_export_hdr.vtu";
        AmrMeshExport::write_vtu(path, &h, 0).expect("write");
        let mut raw = String::new();
        File::open(path).unwrap().read_to_string(&mut raw).unwrap();
        assert!(raw.contains("UnstructuredGrid"));
    }
    #[test]
    fn test_amr_mesh_export_hex_type() {
        let h = make_two_level_hierarchy();
        let path = "/tmp/test_amr_export_hex.vtu";
        AmrMeshExport::write_vtu(path, &h, 0).expect("write");
        let mut raw = String::new();
        File::open(path).unwrap().read_to_string(&mut raw).unwrap();
        assert!(raw.contains("12"));
    }
    #[test]
    fn test_refinement_user_defined() {
        let h = make_two_level_hierarchy();
        let grid = &h.levels[0];
        let crit = AmrRefinementCriteria::new(RefinementMethod::UserDefined {
            field_idx: 0,
            value: 2.0_f64,
        });
        let flags = crit.flag_refinement(grid);
        assert_eq!(flags.len(), 3);
        assert!(!flags[0]);
        assert!(flags[1]);
    }
    #[test]
    fn test_refinement_gradient_based_zero_threshold() {
        let h = make_two_level_hierarchy();
        let grid = &h.levels[1];
        let crit = AmrRefinementCriteria::new(RefinementMethod::GradientBased {
            field_idx: 0,
            threshold: 0.0_f64,
        });
        let flags = crit.flag_refinement(grid);
        assert_eq!(flags.len(), grid.cell_count());
    }
    #[test]
    fn test_refinement_gradient_based_high_threshold() {
        let h = make_two_level_hierarchy();
        let grid = &h.levels[0];
        let crit = AmrRefinementCriteria::new(RefinementMethod::GradientBased {
            field_idx: 0,
            threshold: 1.0e10_f64,
        });
        let flags = crit.flag_refinement(grid);
        assert!(flags.iter().all(|&f| !f));
    }
    #[test]
    fn test_refinement_curvature_high_threshold() {
        let h = make_two_level_hierarchy();
        let grid = &h.levels[0];
        let crit = AmrRefinementCriteria::new(RefinementMethod::CurvatureBased {
            field_idx: 0,
            threshold: 1.0e10_f64,
        });
        let flags = crit.flag_refinement(grid);
        assert!(flags.iter().all(|&f| !f));
    }
    #[test]
    fn test_amr_coarsening_basic() {
        let h = make_two_level_hierarchy();
        let fine = &h.levels[1];
        let coarse_lm = AmrLevel::new(0, 1.0_f64, [0.0; 6]);
        let coarse = AmrCoarsening::coarsen_grid(fine, coarse_lm, 2);
        assert!(!coarse.cells.is_empty());
    }
    #[test]
    fn test_amr_coarsening_averaged_values() {
        let h = make_two_level_hierarchy();
        let fine = &h.levels[1];
        let coarse_lm = AmrLevel::new(0, 1.0_f64, [0.0_f64; 6]);
        let coarse = AmrCoarsening::coarsen_grid(fine, coarse_lm, 2);
        let expected = (1.1_f64 + 1.2_f64 + 1.3_f64 + 1.4_f64) / 4.0_f64;
        let cell = coarse.cells.iter().find(|c| c.coords == [0, 0, 0]).unwrap();
        assert!((cell.data[0] - expected).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_amr_coarsening_field_names_preserved() {
        let h = make_two_level_hierarchy();
        let fine = &h.levels[1];
        let coarse_lm = AmrLevel::new(0, 1.0_f64, [0.0_f64; 6]);
        let coarse = AmrCoarsening::coarsen_grid(fine, coarse_lm, 2);
        assert_eq!(coarse.field_names, vec!["pressure", "density"]);
    }
    #[test]
    fn test_amr_remove_covered_cells() {
        let mut h = AmrHierarchy::new();
        let lm0 = AmrLevel::new(0, 1.0_f64, [0.0; 6]);
        let mut g0 = AmrGrid::new(lm0, vec!["u".into()]);
        g0.add_cell(AmrCell::new(0, [0, 0, 0], vec![1.0_f64]));
        g0.add_cell(AmrCell::new(0, [1, 0, 0], vec![2.0_f64]));
        h.add_level(g0);
        let lm1 = AmrLevel::new(0, 1.0_f64, [0.0; 6]);
        let mut g1 = AmrGrid::new(lm1, vec!["u".into()]);
        g1.add_cell(AmrCell::new(0, [0, 0, 0], vec![1.1_f64]));
        h.add_level(g1);
        let mut coarse = h.levels[0].clone();
        let fine = &h.levels[1];
        let count_before = coarse.cell_count();
        AmrCoarsening::remove_covered_cells(&mut coarse, fine, 1);
        let count_after = coarse.cell_count();
        assert!(count_after < count_before, "covered cell should be removed");
    }
}
