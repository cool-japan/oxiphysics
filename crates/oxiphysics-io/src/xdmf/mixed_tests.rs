// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `Mixed` topology: `nodes_per_element()` is 0, which used to panic through
//! `slice::chunks(0)` in `write_xdmf_unstructured`, `XdmfMeshTimeSeries::to_xml`
//! and `XdmfWriter::write_unstructured_topology`.

use super::mixed::{fixed_rows, mixed_rows};
use super::*;
use std::io::ErrorKind;

/// A triangle (0 1 2) followed by a quadrilateral (0 1 2 3) in XDMF mixed
/// stream layout: `[code, node ids...]` per element.
const TRI_QUAD: [usize; 9] = [4, 0, 1, 2, 5, 0, 1, 2, 3];

fn square_nodes() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
}

fn mixed_step(time: f64, connectivity: Vec<usize>) -> XdmfMeshStep {
    XdmfMeshStep {
        time,
        nodes: square_nodes(),
        connectivity,
        topology: XdmfTopologyType::Mixed,
        node_scalars: vec![],
        node_vectors: vec![],
    }
}

// ── stream walker ────────────────────────────────────────────────────────

#[test]
fn counts_a_triangle_and_a_quadrilateral() {
    assert_eq!(xdmf_mixed_element_count(&TRI_QUAD).unwrap(), 2);
}

#[test]
fn empty_stream_has_no_elements() {
    assert_eq!(xdmf_mixed_element_count(&[]).unwrap(), 0);
    assert!(mixed_rows(&[]).unwrap().is_empty());
}

#[test]
fn counts_every_fixed_size_cell_type() {
    // Triangle 3, Quadrilateral 4, Tetrahedron 4, Pyramid 5, Wedge 6, Hexahedron 8.
    let stream = [
        4, 0, 1, 2, //
        5, 0, 1, 2, 3, //
        6, 0, 1, 2, 3, //
        7, 0, 1, 2, 3, 4, //
        8, 0, 1, 2, 3, 4, 5, //
        9, 0, 1, 2, 3, 4, 5, 6, 7,
    ];
    assert_eq!(xdmf_mixed_element_count(&stream).unwrap(), 6);
    let rows = mixed_rows(&stream).unwrap();
    assert_eq!(rows[0], "4 0 1 2");
    assert_eq!(rows[3], "7 0 1 2 3 4");
    assert_eq!(rows[5], "9 0 1 2 3 4 5 6 7");
}

#[test]
fn polytypes_carry_an_explicit_node_count() {
    // Polyvertex (2 nodes), Polyline (3 nodes), Polygon (5 nodes).
    let stream = [
        1, 2, 10, 11, //
        2, 3, 0, 1, 2, //
        3, 5, 0, 1, 2, 3, 4,
    ];
    assert_eq!(xdmf_mixed_element_count(&stream).unwrap(), 3);
    let rows = mixed_rows(&stream).unwrap();
    assert_eq!(rows, vec!["1 2 10 11", "2 3 0 1 2", "3 5 0 1 2 3 4"]);
}

#[test]
fn rejects_unknown_cell_type_codes() {
    for code in [0usize, 10, 34, 112] {
        let err = xdmf_mixed_element_count(&[4, 0, 1, 2, code, 0, 1]).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidInput, "code {code}");
        assert!(err.to_string().contains(&format!("code {code}")), "{err}");
        assert!(err.to_string().contains("offset 4"), "{err}");
    }
}

#[test]
fn rejects_an_element_that_runs_past_the_end() {
    // The quadrilateral needs four node ids but only three follow.
    let err = xdmf_mixed_element_count(&[4, 0, 1, 2, 5, 0, 1, 2]).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidInput);
    assert!(err.to_string().contains("runs past the end"), "{err}");
    // A lone cell-type code with no nodes at all.
    assert!(xdmf_mixed_element_count(&[9]).is_err());
    // A polytype without its node count.
    assert!(xdmf_mixed_element_count(&[2]).is_err());
    // A polytype whose count exceeds the remaining entries (and cannot overflow).
    assert!(xdmf_mixed_element_count(&[3, usize::MAX, 0, 1]).is_err());
}

#[test]
fn rejects_polytypes_with_too_few_nodes() {
    assert!(xdmf_mixed_element_count(&[1, 0]).is_err());
    assert!(xdmf_mixed_element_count(&[2, 1, 7]).is_err());
    assert!(xdmf_mixed_element_count(&[3, 2, 7, 8]).is_err());
}

#[test]
fn fixed_rows_group_by_node_count() {
    assert_eq!(fixed_rows(&[0, 1, 2, 3, 4, 5], 3), vec!["0 1 2", "3 4 5"]);
    // A ragged tail is kept as a shorter last row.
    assert_eq!(fixed_rows(&[0, 1, 2, 3], 3), vec!["0 1 2", "3"]);
    // A zero node count must not panic.
    assert_eq!(fixed_rows(&[0, 1], 0), vec!["0", "1"]);
}

// ── write_xdmf_unstructured ──────────────────────────────────────────────

#[test]
fn write_unstructured_mixed_tri_and_quad() {
    let mut buf = Vec::new();
    write_xdmf_unstructured(
        &mut buf,
        XdmfTopologyType::Mixed,
        &square_nodes(),
        &TRI_QUAD,
        &[],
    )
    .unwrap();
    let s = String::from_utf8(buf).unwrap();
    assert!(
        s.contains("<Topology TopologyType=\"Mixed\" NumberOfElements=\"2\">"),
        "{s}"
    );
    assert!(
        s.contains("<DataItem Format=\"XML\" DataType=\"Int\" Dimensions=\"9\">"),
        "{s}"
    );
    assert!(s.contains("          4 0 1 2\n"), "{s}");
    assert!(s.contains("          5 0 1 2 3\n"), "{s}");
    assert!(
        s.contains("<DataItem Format=\"XML\" Dimensions=\"4 3\">"),
        "{s}"
    );
}

#[test]
fn write_unstructured_mixed_rejects_a_malformed_stream_before_writing() {
    let mut buf = Vec::new();
    let err = write_xdmf_unstructured(
        &mut buf,
        XdmfTopologyType::Mixed,
        &square_nodes(),
        &[4, 0, 1, 2, 77, 0],
        &[],
    )
    .unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidInput);
    assert!(
        err.to_string().contains("unknown cell type code 77"),
        "{err}"
    );
    assert!(
        buf.is_empty(),
        "nothing may be written for an invalid stream"
    );
}

#[test]
fn write_unstructured_fixed_topologies_are_unchanged() {
    let mut buf = Vec::new();
    write_xdmf_unstructured(
        &mut buf,
        XdmfTopologyType::Triangle,
        &square_nodes(),
        &[0, 1, 2, 0, 2, 3],
        &[],
    )
    .unwrap();
    let s = String::from_utf8(buf).unwrap();
    assert!(s.contains("NumberOfElements=\"2\""), "{s}");
    assert!(
        s.contains("        <DataItem Format=\"XML\" Dimensions=\"2 3\">\n"),
        "{s}"
    );
    assert!(!s.contains("DataType"), "{s}");
    assert!(s.contains("          0 1 2\n          0 2 3\n"), "{s}");
}

// ── XdmfMeshTimeSeries / XdmfMeshStep ────────────────────────────────────

#[test]
fn mesh_step_counts_mixed_elements_by_walking_the_stream() {
    assert_eq!(mixed_step(0.0, TRI_QUAD.to_vec()).n_elements(), 2);
    // A malformed stream counts as no elements instead of panicking.
    assert_eq!(mixed_step(0.0, vec![4, 0, 1]).n_elements(), 0);
    assert_eq!(mixed_step(0.0, vec![]).n_elements(), 0);
}

#[test]
fn mesh_series_to_xml_writes_a_mixed_step() {
    let mut ms = XdmfMeshTimeSeries::new();
    ms.add_step(mixed_step(0.5, TRI_QUAD.to_vec()));
    let xml = ms.to_xml().unwrap();
    assert!(
        xml.contains("<Topology TopologyType=\"Mixed\" NumberOfElements=\"2\">"),
        "{xml}"
    );
    assert!(
        xml.contains("<DataItem Format=\"XML\" DataType=\"Int\" Dimensions=\"9\">"),
        "{xml}"
    );
    assert!(xml.contains("            4 0 1 2\n"), "{xml}");
    assert!(xml.contains("            5 0 1 2 3\n"), "{xml}");
    assert!(xml.contains("<Time Value=\"0.5\"/>"), "{xml}");
}

#[test]
fn mesh_series_to_xml_reports_the_malformed_step() {
    let mut ms = XdmfMeshTimeSeries::new();
    ms.add_step(mixed_step(0.0, TRI_QUAD.to_vec()));
    ms.add_step(mixed_step(1.5, vec![4, 0, 1, 2, 6, 0, 1]));
    let err = ms.to_xml().unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidInput);
    assert!(err.to_string().contains("mesh step 1"), "{err}");
    assert!(err.to_string().contains("t = 1.5"), "{err}");
}

#[test]
fn mesh_series_with_fixed_and_mixed_steps() {
    let mut ms = XdmfMeshTimeSeries::new();
    ms.add_step(XdmfMeshStep {
        time: 0.0,
        nodes: square_nodes(),
        connectivity: vec![0, 1, 2],
        topology: XdmfTopologyType::Triangle,
        node_scalars: vec![],
        node_vectors: vec![],
    });
    ms.add_step(mixed_step(1.0, TRI_QUAD.to_vec()));
    let xml = ms.to_xml().unwrap();
    assert!(xml.contains("<Topology TopologyType=\"Triangle\" NumberOfElements=\"1\">"));
    assert!(xml.contains("<Topology TopologyType=\"Mixed\" NumberOfElements=\"2\">"));
    // The mixed step has more elements, so it is the peak step.
    assert_eq!(peak_element_step(&ms), Some(1));
}

// ── XdmfWriter::write_unstructured_topology ──────────────────────────────

#[test]
fn writer_mixed_topology_is_written_flat() {
    let mut w = XdmfWriter::new();
    w.write_unstructured_topology("Mixed", 2, &TRI_QUAD, 0);
    let xml = w.finish();
    assert!(
        xml.contains("<Topology TopologyType=\"Mixed\" NumberOfElements=\"2\">"),
        "{xml}"
    );
    assert!(
        xml.contains("<DataItem Format=\"XML\" DataType=\"Int\" Dimensions=\"9\">"),
        "{xml}"
    );
    assert!(xml.contains("      4 0 1 2\n"), "{xml}");
    assert!(xml.contains("      5 0 1 2 3\n"), "{xml}");
    assert!(xml.contains("</DataItem>\n  </Topology>"), "{xml}");
}

#[test]
fn writer_mixed_topology_does_not_panic_on_a_malformed_stream() {
    let mut w = XdmfWriter::new();
    w.write_unstructured_topology("Mixed", 1, &[4, 0, 1], 0);
    let xml = w.finish();
    // The stream cannot be split into elements, so it is written as one line.
    assert!(xml.contains("Dimensions=\"3\""), "{xml}");
    assert!(xml.contains("4 0 1\n"), "{xml}");
    assert!(xml.contains("</Topology>"), "{xml}");
}

#[test]
fn writer_fixed_topology_is_unchanged() {
    let mut w = XdmfWriter::new();
    w.write_unstructured_topology("Triangle", 2, &[0, 1, 2, 0, 2, 3], 3);
    let xml = w.finish();
    assert!(xml.contains("Dimensions=\"2 3\""), "{xml}");
    assert!(xml.contains("0 1 2\n"), "{xml}");
    assert!(xml.contains("0 2 3\n"), "{xml}");
    assert!(!xml.contains("DataType"), "{xml}");
}
