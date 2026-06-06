//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

// Types are re-exported via mod.rs; import in test modules via use super::*.

/// Count how many `<Grid ...>` elements appear in an XDMF XML string.
pub fn count_xdmf_grids(xml: &str) -> usize {
    xml.matches("<Grid ").count()
}
/// Extract all `Name="..."` attribute values from XDMF XML.
pub fn extract_xdmf_names(xml: &str) -> Vec<String> {
    let mut names = Vec::new();
    let mut rest = xml;
    while let Some(pos) = rest.find("Name=\"") {
        rest = &rest[pos + 6..];
        if let Some(end) = rest.find('"') {
            names.push(rest[..end].to_string());
            rest = &rest[end + 1..];
        } else {
            break;
        }
    }
    names
}
#[cfg(test)]
mod tests_xdmf_additions {
    use crate::xdmf::*;
    #[test]
    fn field_scalar_entry_count() {
        let f = XdmfFieldDescriptor::scalar("p", vec![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(f.entry_count(), 4);
    }
    #[test]
    fn field_vector_entry_count() {
        let data: Vec<f64> = (0..12).map(|x| x as f64).collect();
        let f = XdmfFieldDescriptor::vector("vel", data);
        assert_eq!(f.entry_count(), 4);
    }
    #[test]
    fn field_max_abs_positive() {
        let f = XdmfFieldDescriptor::scalar("p", vec![-5.0, 2.0, 3.0]);
        assert!((f.max_abs() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn field_min_value_correct() {
        let f = XdmfFieldDescriptor::scalar("p", vec![3.0, -1.5, 2.0]);
        assert!((f.min_value() - (-1.5)).abs() < 1e-12);
    }
    #[test]
    fn field_max_value_correct() {
        let f = XdmfFieldDescriptor::scalar("p", vec![3.0, -1.5, 7.0]);
        assert!((f.max_value() - 7.0).abs() < 1e-12);
    }
    #[test]
    fn field_l2_norm_basic() {
        let f = XdmfFieldDescriptor::scalar("p", vec![3.0, 4.0]);
        let norm = f.data_lp_norm(2.0);
        assert!((norm - 5.0).abs() < 1e-10);
    }
    #[test]
    fn field_zero_components_returns_zero_entries() {
        let f = XdmfFieldDescriptor {
            name: "x".to_string(),
            attribute_type: "Scalar".to_string(),
            center: "Node".to_string(),
            data: vec![1.0, 2.0],
            n_components: 0,
        };
        assert_eq!(f.entry_count(), 0);
    }
    #[test]
    fn write_timestep_fields_produces_grid_tag() {
        let nodes = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let fields = vec![XdmfFieldDescriptor::scalar("pressure", vec![1.0, 2.0, 3.0])];
        let mut buf = Vec::new();
        write_xdmf_timestep_fields(&mut buf, 0.5, 3, "Triangle", 1, &nodes, &fields).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("<Grid"));
        assert!(s.contains("pressure"));
        assert!(s.contains("0.5"));
    }
    #[test]
    fn write_timestep_fields_no_fields_still_valid() {
        let nodes = vec![[0.0_f64; 3]; 2];
        let mut buf = Vec::new();
        write_xdmf_timestep_fields(&mut buf, 1.0, 2, "Polyline", 1, &nodes, &[]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Geometry"));
    }
    #[test]
    fn patch_element_count() {
        let p = XdmfMeshPatch::new("inlet", vec![0, 1, 2, 3]);
        assert_eq!(p.element_count(), 4);
    }
    #[test]
    fn patch_contains_element_true() {
        let p = XdmfMeshPatch::new("wall", vec![5, 10, 15]);
        assert!(p.contains_element(10));
    }
    #[test]
    fn patch_contains_element_false() {
        let p = XdmfMeshPatch::new("wall", vec![5, 10, 15]);
        assert!(!p.contains_element(7));
    }
    #[test]
    fn patch_merge_no_duplicates() {
        let mut p1 = XdmfMeshPatch::new("a", vec![1, 2, 3]);
        let p2 = XdmfMeshPatch::new("b", vec![2, 4, 5]);
        p1.merge(&p2);
        assert_eq!(p1.element_count(), 5);
    }
    #[test]
    fn patch_to_debug_string_contains_name() {
        let p = XdmfMeshPatch::new("outlet", vec![0]);
        let s = p.to_debug_string();
        assert!(s.contains("outlet"));
    }
    #[test]
    fn patch_element_map_first_match_wins() {
        let patches = vec![
            XdmfMeshPatch::new("a", vec![0, 1]),
            XdmfMeshPatch::new("b", vec![1, 2]),
        ];
        let map = patch_element_map(&patches);
        assert_eq!(map.get(&1).map(|s| s.as_str()), Some("a"));
        assert_eq!(map.get(&2).map(|s| s.as_str()), Some("b"));
    }
    #[test]
    fn grid_entry_bounding_box_unit_triangle() {
        let entry = XdmfGridEntry {
            name: "t".to_string(),
            time: 0.0,
            nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            connectivity: vec![0, 1, 2],
            topology_type: "Triangle".to_string(),
        };
        let (lo, hi) = entry.bounding_box();
        assert!((lo[0] - 0.0).abs() < 1e-12);
        assert!((hi[0] - 1.0).abs() < 1e-12);
        assert!((hi[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn grid_entry_centroid_equilateral() {
        let entry = XdmfGridEntry {
            name: "t".to_string(),
            time: 0.0,
            nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            connectivity: vec![0, 1, 2],
            topology_type: "Triangle".to_string(),
        };
        let c = entry.centroid();
        assert!((c[0] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn grid_entry_empty_nodes_centroid_zero() {
        let entry = XdmfGridEntry {
            name: "e".to_string(),
            time: 0.0,
            nodes: vec![],
            connectivity: vec![],
            topology_type: "Triangle".to_string(),
        };
        let c = entry.centroid();
        assert_eq!(c, [0.0; 3]);
    }
    #[test]
    fn domain_collection_total_node_count() {
        let mut dc = XdmfDomainCollection::new();
        dc.add_grid(XdmfGridEntry {
            name: "a".to_string(),
            time: 0.0,
            nodes: vec![[0.0; 3]; 4],
            connectivity: vec![],
            topology_type: "Triangle".to_string(),
        });
        dc.add_grid(XdmfGridEntry {
            name: "b".to_string(),
            time: 1.0,
            nodes: vec![[0.0; 3]; 6],
            connectivity: vec![],
            topology_type: "Triangle".to_string(),
        });
        assert_eq!(dc.total_node_count(), 10);
    }
    #[test]
    fn domain_collection_find_grid_hit() {
        let mut dc = XdmfDomainCollection::new();
        dc.add_grid(XdmfGridEntry {
            name: "target".to_string(),
            time: 0.5,
            nodes: vec![],
            connectivity: vec![],
            topology_type: "Polyvertex".to_string(),
        });
        assert!(dc.find_grid("target").is_some());
    }
    #[test]
    fn domain_collection_find_grid_miss() {
        let dc = XdmfDomainCollection::new();
        assert!(dc.find_grid("nonexistent").is_none());
    }
    #[test]
    fn domain_collection_time_range_empty() {
        let dc = XdmfDomainCollection::new();
        assert_eq!(dc.time_range(), (0.0, 0.0));
    }
    #[test]
    fn domain_collection_time_range_multiple_grids() {
        let mut dc = XdmfDomainCollection::new();
        for t in [0.5_f64, 1.0, 2.5] {
            dc.add_grid(XdmfGridEntry {
                name: format!("g{t}"),
                time: t,
                nodes: vec![],
                connectivity: vec![],
                topology_type: "Triangle".to_string(),
            });
        }
        let (t_min, t_max) = dc.time_range();
        assert!((t_min - 0.5).abs() < 1e-12);
        assert!((t_max - 2.5).abs() < 1e-12);
    }
    #[test]
    fn domain_collection_write_xml_valid_structure() {
        let mut dc = XdmfDomainCollection::new();
        dc.add_grid(XdmfGridEntry {
            name: "step0".to_string(),
            time: 0.0,
            nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            connectivity: vec![0, 1, 2],
            topology_type: "Triangle".to_string(),
        });
        let mut buf = Vec::new();
        dc.write_xml(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("<?xml"));
        assert!(s.contains("</Xdmf>"));
        assert!(s.contains("step0"));
    }
    #[test]
    fn format_geometry_inline_two_nodes() {
        let nodes = vec![[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let s = format_xdmf_geometry_inline(&nodes);
        assert!(s.contains("1 2 3") || s.contains("1.0 2.0 3.0") || s.contains("1 2"));
    }
    #[test]
    fn format_data_inline_basic() {
        let s = format_xdmf_data_inline(&[1.0, 2.0, 3.0]);
        assert_eq!(s, "1 2 3");
    }
    #[test]
    fn xdmf_hdf5_dataitem_contains_path() {
        let s = xdmf_hdf5_dataitem("sim.h5", "/data/pressure", "100", "Float");
        assert!(s.contains("sim.h5:/data/pressure"));
        assert!(s.contains("HDF"));
    }
    #[test]
    fn validate_xdmf_structure_valid() {
        let xml = "<?xml version=\"1.0\"?><Xdmf><Domain></Domain></Xdmf>";
        assert!(validate_xdmf_structure(xml).is_ok());
    }
    #[test]
    fn validate_xdmf_structure_missing_domain() {
        let xml = "<?xml version=\"1.0\"?><Xdmf></Xdmf>";
        assert!(validate_xdmf_structure(xml).is_err());
    }
    #[test]
    fn indent_xdmf_adds_spaces() {
        let s = indent_xdmf("<foo/>", 2);
        assert!(s.starts_with("    "));
    }
    #[test]
    fn count_xdmf_grids_two_grids() {
        let xml = "<Grid Name=\"a\"><Grid Name=\"b\"></Grid></Grid>";
        assert_eq!(count_xdmf_grids(xml), 2);
    }
    #[test]
    fn extract_xdmf_names_finds_all() {
        let xml = "<Grid Name=\"foo\"><Attribute Name=\"bar\"/></Grid>";
        let names = extract_xdmf_names(xml);
        assert!(names.contains(&"foo".to_string()));
        assert!(names.contains(&"bar".to_string()));
    }
}
