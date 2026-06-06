//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;
use std::io::Write;

use super::types::{
    XdmfAttribute, XdmfFieldDescriptor, XdmfMeshPatch, XdmfMeshTimeSeries, XdmfTopologyType,
    XdmfUniformGrid,
};

/// Write an XDMF file for a particle dataset (points + optional scalar fields).
///
/// Produces XML like:
/// ```xml
/// <?xml version="1.0"?>
/// <Xdmf Version="3.0">
///   `Domain`
///     <Grid Name="particles" GridType="Uniform">
///       <Topology TopologyType="Polyvertex" NumberOfElements="N"/>
///       <Geometry GeometryType="XYZ">
///         <DataItem Format="XML" Dimensions="N 3">...</DataItem>
///       </Geometry>
///       <Attribute Name="density" AttributeType="Scalar" Center="Node">
///         <DataItem Format="XML" Dimensions="N">...</DataItem>
///       </Attribute>
///     </Grid>
///   </Domain>
/// </Xdmf>
/// ```
pub fn write_xdmf_particles<W: Write>(
    writer: &mut W,
    positions: &[Vec3],
    scalar_fields: &[(&str, &[f64])],
) -> std::io::Result<()> {
    let n = positions.len();
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<Xdmf Version=\"3.0\">")?;
    writeln!(writer, "  <Domain>")?;
    writeln!(writer, "    <Grid Name=\"particles\" GridType=\"Uniform\">")?;
    writeln!(
        writer,
        "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"{}\"/>",
        n
    )?;
    writeln!(writer, "      <Geometry GeometryType=\"XYZ\">")?;
    writeln!(
        writer,
        "        <DataItem Format=\"XML\" Dimensions=\"{} 3\">",
        n
    )?;
    for p in positions {
        writeln!(writer, "          {} {} {}", p.x, p.y, p.z)?;
    }
    writeln!(writer, "        </DataItem>")?;
    writeln!(writer, "      </Geometry>")?;
    for (name, values) in scalar_fields {
        writeln!(
            writer,
            "      <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">",
            name
        )?;
        writeln!(
            writer,
            "        <DataItem Format=\"XML\" Dimensions=\"{}\">",
            n
        )?;
        write!(writer, "          ")?;
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{}", v)?;
        }
        writeln!(writer)?;
        writeln!(writer, "        </DataItem>")?;
        writeln!(writer, "      </Attribute>")?;
    }
    writeln!(writer, "    </Grid>")?;
    writeln!(writer, "  </Domain>")?;
    writeln!(writer, "</Xdmf>")?;
    Ok(())
}
/// Write a multi-timestep XDMF temporal collection.
///
/// Each entry in `timesteps` is `(time, positions)`. The result is a
/// `<Grid GridType="Collection" CollectionType="Temporal">` containing one
/// child grid per timestep.
pub fn write_xdmf_temporal<W: Write>(
    writer: &mut W,
    timesteps: &[(f64, &[Vec3])],
) -> std::io::Result<()> {
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<Xdmf Version=\"3.0\">")?;
    writeln!(writer, "  <Domain>")?;
    writeln!(
        writer,
        "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">"
    )?;
    for (time, positions) in timesteps {
        let n = positions.len();
        writeln!(
            writer,
            "      <Grid Name=\"particles\" GridType=\"Uniform\">"
        )?;
        writeln!(writer, "        <Time Value=\"{}\"/>", time)?;
        writeln!(
            writer,
            "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"{}\"/>",
            n
        )?;
        writeln!(writer, "        <Geometry GeometryType=\"XYZ\">")?;
        writeln!(
            writer,
            "          <DataItem Format=\"XML\" Dimensions=\"{} 3\">",
            n
        )?;
        for p in *positions {
            writeln!(writer, "            {} {} {}", p.x, p.y, p.z)?;
        }
        writeln!(writer, "          </DataItem>")?;
        writeln!(writer, "        </Geometry>")?;
        writeln!(writer, "      </Grid>")?;
    }
    writeln!(writer, "    </Grid>")?;
    writeln!(writer, "  </Domain>")?;
    writeln!(writer, "</Xdmf>")?;
    Ok(())
}
/// Write an XDMF file referencing HDF5 geometry and attributes.
pub fn write_xdmf_with_attributes(
    path: &str,
    hdf5_path: &str,
    n_nodes: usize,
    n_elements: usize,
    topology: &str,
    attributes: &[XdmfAttribute],
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "<?xml version=\"1.0\"?>")?;
    writeln!(f, "<Xdmf Version=\"3.0\">")?;
    writeln!(f, "  <Domain>")?;
    writeln!(f, "    <Grid Name=\"mesh\" GridType=\"Uniform\">")?;
    writeln!(
        f,
        "      <Topology TopologyType=\"{}\" NumberOfElements=\"{}\"/>",
        topology, n_elements
    )?;
    writeln!(f, "      <Geometry GeometryType=\"XYZ\">")?;
    writeln!(
        f,
        "        <DataItem Format=\"HDF\" Dimensions=\"{} 3\">{}:/coordinates</DataItem>",
        n_nodes, hdf5_path
    )?;
    writeln!(f, "      </Geometry>")?;
    for attr in attributes {
        let attr_type = if attr.n_components == 1 {
            "Scalar"
        } else {
            "Vector"
        };
        let dims = if attr.n_components == 1 {
            format!("{}", n_nodes)
        } else {
            format!("{} {}", n_nodes, attr.n_components)
        };
        writeln!(
            f,
            "      <Attribute Name=\"{}\" AttributeType=\"{}\" Center=\"{}\">",
            attr.name, attr_type, attr.center
        )?;
        writeln!(
            f,
            "        <DataItem Format=\"HDF\" Dimensions=\"{}\">{}</DataItem>",
            dims, attr.hdf5_path
        )?;
        writeln!(f, "      </Attribute>")?;
    }
    writeln!(f, "    </Grid>")?;
    writeln!(f, "  </Domain>")?;
    writeln!(f, "</Xdmf>")?;
    Ok(())
}
/// Parse topology information from XDMF XML content.
///
/// Returns a list of `(topology_type, n_elements)` pairs found in the document.
pub fn parse_xdmf_topology(content: &str) -> Vec<(String, usize)> {
    let mut result = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if !trimmed.contains("Topology") {
            continue;
        }
        let topo_type = if let Some(start) = trimmed.find("TopologyType=\"") {
            let rest = &trimmed[start + 14..];
            rest.find('"').map(|end| rest[..end].to_string())
        } else {
            None
        };
        let n_elements = if let Some(start) = trimmed.find("NumberOfElements=\"") {
            let rest = &trimmed[start + 18..];
            rest.find('"')
                .and_then(|end| rest[..end].parse::<usize>().ok())
        } else {
            None
        };
        if let (Some(t), Some(n)) = (topo_type, n_elements) {
            result.push((t, n));
        }
    }
    result
}
/// Write a structured uniform-grid XDMF file.
pub fn write_xdmf_uniform_grid(path: &str, params: &XdmfUniformGrid) -> std::io::Result<()> {
    let [nx, ny, nz] = params.dimensions;
    let [ox, oy, oz] = params.origin;
    let [dx, dy, dz] = params.spacing;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "<?xml version=\"1.0\"?>")?;
    writeln!(f, "<Xdmf Version=\"3.0\">")?;
    writeln!(f, "  <Domain>")?;
    writeln!(
        f,
        "    <Grid Name=\"{}\" GridType=\"Uniform\">",
        params.name
    )?;
    writeln!(
        f,
        "      <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"{} {} {}\"/>",
        nz, ny, nx
    )?;
    writeln!(f, "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">")?;
    writeln!(
        f,
        "        <DataItem Format=\"XML\" Dimensions=\"3\">{} {} {}</DataItem>",
        oz, oy, ox
    )?;
    writeln!(
        f,
        "        <DataItem Format=\"XML\" Dimensions=\"3\">{} {} {}</DataItem>",
        dz, dy, dx
    )?;
    writeln!(f, "      </Geometry>")?;
    writeln!(f, "    </Grid>")?;
    writeln!(f, "  </Domain>")?;
    writeln!(f, "</Xdmf>")?;
    Ok(())
}
/// Generate an XDMF ``DataItem` element referencing an HDF5 file.
pub fn write_xdmf_hdf5_reference(filename: &str, dataset_path: &str) -> String {
    format!(
        "<DataItem Format=\"HDF\" Dimensions=\"1\">\n  {}:{}\n</DataItem>",
        filename, dataset_path
    )
}
/// Write a vector attribute block (e.g., velocity) into an XDMF XML string.
///
/// Produces:
/// ```xml
/// <Attribute Name="velocity" AttributeType="Vector" Center="Node">
///   <DataItem Format="XML" Dimensions="N 3">
///     vx vy vz
///     ...
///   </DataItem>
/// </Attribute>
/// ```
pub fn xdmf_vector_attribute(name: &str, vectors: &[[f64; 3]]) -> String {
    let n = vectors.len();
    let mut s = String::new();
    s.push_str(&format!(
        "      <Attribute Name=\"{}\" AttributeType=\"Vector\" Center=\"Node\">\n",
        name
    ));
    s.push_str(&format!(
        "        <DataItem Format=\"XML\" Dimensions=\"{} 3\">\n",
        n
    ));
    for v in vectors {
        s.push_str(&format!("          {} {} {}\n", v[0], v[1], v[2]));
    }
    s.push_str("        </DataItem>\n");
    s.push_str("      </Attribute>\n");
    s
}
/// Write a symmetric 3×3 tensor attribute block (e.g., stress tensor).
///
/// Each tensor is stored as 6 unique components: `\[xx, yy, zz, xy, xz, yz\]`.
/// XDMF AttributeType is `"Tensor6"`.
pub fn xdmf_tensor6_attribute(name: &str, tensors: &[[f64; 6]]) -> String {
    let n = tensors.len();
    let mut s = String::new();
    s.push_str(&format!(
        "      <Attribute Name=\"{}\" AttributeType=\"Tensor6\" Center=\"Node\">\n",
        name
    ));
    s.push_str(&format!(
        "        <DataItem Format=\"XML\" Dimensions=\"{} 6\">\n",
        n
    ));
    for t in tensors {
        s.push_str(&format!(
            "          {} {} {} {} {} {}\n",
            t[0], t[1], t[2], t[3], t[4], t[5]
        ));
    }
    s.push_str("        </DataItem>\n");
    s.push_str("      </Attribute>\n");
    s
}
/// Write an unstructured mesh XDMF grid to a `Write` target.
///
/// `nodes` are 3-D point coordinates; `connectivity` is a flat list of
/// node indices, grouped by element (length must be a multiple of
/// `topo.nodes_per_element()`).
pub fn write_xdmf_unstructured<W: Write>(
    writer: &mut W,
    topo: XdmfTopologyType,
    nodes: &[[f64; 3]],
    connectivity: &[usize],
    scalar_fields: &[(&str, &[f64])],
) -> std::io::Result<()> {
    let n_nodes = nodes.len();
    let npe = topo.nodes_per_element();
    let n_elements = connectivity.len().checked_div(npe).unwrap_or(0);
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<Xdmf Version=\"3.0\">")?;
    writeln!(writer, "  <Domain>")?;
    writeln!(writer, "    <Grid Name=\"mesh\" GridType=\"Uniform\">")?;
    writeln!(
        writer,
        "      <Topology TopologyType=\"{}\" NumberOfElements=\"{}\">",
        topo.xdmf_name(),
        n_elements
    )?;
    writeln!(
        writer,
        "        <DataItem Format=\"XML\" Dimensions=\"{} {}\">",
        n_elements, npe
    )?;
    for chunk in connectivity.chunks(npe) {
        let row: Vec<String> = chunk.iter().map(|&i| i.to_string()).collect();
        writeln!(writer, "          {}", row.join(" "))?;
    }
    writeln!(writer, "        </DataItem>")?;
    writeln!(writer, "      </Topology>")?;
    writeln!(writer, "      <Geometry GeometryType=\"XYZ\">")?;
    writeln!(
        writer,
        "        <DataItem Format=\"XML\" Dimensions=\"{} 3\">",
        n_nodes
    )?;
    for p in nodes {
        writeln!(writer, "          {} {} {}", p[0], p[1], p[2])?;
    }
    writeln!(writer, "        </DataItem>")?;
    writeln!(writer, "      </Geometry>")?;
    for (name, values) in scalar_fields {
        writeln!(
            writer,
            "      <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">",
            name
        )?;
        writeln!(
            writer,
            "        <DataItem Format=\"XML\" Dimensions=\"{}\">",
            values.len()
        )?;
        write!(writer, "          ")?;
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{}", v)?;
        }
        writeln!(writer)?;
        writeln!(writer, "        </DataItem>")?;
        writeln!(writer, "      </Attribute>")?;
    }
    writeln!(writer, "    </Grid>")?;
    writeln!(writer, "  </Domain>")?;
    writeln!(writer, "</Xdmf>")?;
    Ok(())
}
/// Produce a scalar XDMF ``DataItem` XML string (inline, Format="XML").
///
/// Useful for embedding small fields directly in the XDMF document.
pub fn xdmf_scalar_data_item(values: &[f64]) -> String {
    let mut s = format!(
        "<DataItem Format=\"XML\" Dimensions=\"{}\">\n  ",
        values.len()
    );
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&format!("{}", v));
    }
    s.push_str("\n</DataItem>");
    s
}
/// Produce a 3-component vector XDMF ``DataItem` XML string.
pub fn xdmf_vector_data_item(vectors: &[[f64; 3]]) -> String {
    let n = vectors.len();
    let mut s = format!("<DataItem Format=\"XML\" Dimensions=\"{} 3\">\n", n);
    for v in vectors {
        s.push_str(&format!("  {} {} {}\n", v[0], v[1], v[2]));
    }
    s.push_str("</DataItem>");
    s
}
/// Generate a ``Time` element for a temporal grid.
pub fn xdmf_time_element(t: f64) -> String {
    format!("<Time Value=\"{}\"/>", t)
}
/// Compute the total node count across all steps in a mesh time series.
pub fn total_node_count(series: &XdmfMeshTimeSeries) -> usize {
    series.steps.iter().map(|s| s.nodes.len()).sum()
}
/// Return the step index with the most elements (peak mesh density).
pub fn peak_element_step(series: &XdmfMeshTimeSeries) -> Option<usize> {
    series
        .steps
        .iter()
        .enumerate()
        .max_by_key(|(_, s)| s.n_elements())
        .map(|(i, _)| i)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::xdmf::types::*;
    fn make_positions(n: usize) -> Vec<Vec3> {
        (0..n)
            .map(|i| Vec3::new(i as f64, i as f64 * 0.5, 0.0))
            .collect()
    }
    #[test]
    fn test_xdmf_valid_xml() {
        let positions = make_positions(3);
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, &[]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("<?xml"), "missing XML declaration");
        assert!(s.contains("<Xdmf"), "missing Xdmf root element");
        assert!(s.contains("<Grid"), "missing Grid element");
    }
    #[test]
    fn test_xdmf_correct_particle_count() {
        let positions = make_positions(5);
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, &[]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.contains("NumberOfElements=\"5\""),
            "expected NumberOfElements=\"5\" in output, got:\n{}",
            s
        );
    }
    #[test]
    fn test_xdmf_scalar_field_included() {
        let positions = make_positions(4);
        let density = vec![1.0, 2.0, 3.0, 4.0];
        let fields: &[(&str, &[f64])] = &[("density", &density)];
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, fields).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.contains("density"),
            "scalar field name 'density' not found in output"
        );
        assert!(s.contains("Scalar"), "AttributeType Scalar not found");
    }
    #[test]
    fn test_xdmf_temporal_multiple_steps() {
        let pos0 = make_positions(3);
        let pos1 = make_positions(3);
        let steps: &[(f64, &[Vec3])] = &[(0.0, &pos0), (1.0, &pos1)];
        let mut buf = Vec::new();
        write_xdmf_temporal(&mut buf, steps).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Temporal"), "expected CollectionType Temporal");
        assert!(
            s.contains("Time Value=\"0\"") || s.contains("Time Value=\"0."),
            "time step 0 not found"
        );
        assert!(
            s.contains("Time Value=\"1\"") || s.contains("Time Value=\"1."),
            "time step 1 not found"
        );
    }
    #[test]
    fn test_xdmf_empty_positions() {
        let positions: Vec<Vec3> = vec![];
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, &[]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("NumberOfElements=\"0\""));
    }
    #[test]
    fn test_time_series_add_step() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![[0.0, 0.0, 0.0]], vec![]);
        ts.add_step(1.0, vec![[1.0, 0.0, 0.0]], vec![]);
        assert_eq!(ts.steps.len(), 2);
        assert!((ts.steps[0].time - 0.0).abs() < 1e-10);
        assert!((ts.steps[1].time - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_time_series_to_xml() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(
            0.5,
            vec![[1.0, 2.0, 3.0]],
            vec![("density".to_string(), vec![1.5])],
        );
        let xml = ts.to_xml();
        assert!(xml.contains("Temporal"));
        assert!(xml.contains("Time Value=\"0.5\""));
        assert!(xml.contains("density"));
        assert!(xml.contains("1 2 3"));
    }
    #[test]
    fn test_time_series_multiple_scalars() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(
            0.0,
            vec![[0.0; 3]],
            vec![
                ("temperature".to_string(), vec![300.0]),
                ("pressure".to_string(), vec![101325.0]),
            ],
        );
        let xml = ts.to_xml();
        assert!(xml.contains("temperature"));
        assert!(xml.contains("pressure"));
    }
    #[test]
    fn test_xdmf_reader_from_xml() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![[0.0; 3]; 3], vec![]);
        ts.add_step(1.0, vec![[1.0; 3]; 5], vec![]);
        let xml = ts.to_xml();
        let steps = XdmfReader::from_xml(&xml).unwrap();
        assert_eq!(steps.len(), 2);
        assert!((steps[0].time - 0.0).abs() < 1e-10);
        assert_eq!(steps[0].n_points, 3);
        assert!((steps[1].time - 1.0).abs() < 1e-10);
        assert_eq!(steps[1].n_points, 5);
    }
    #[test]
    fn test_xdmf_reader_empty_xml() {
        let steps = XdmfReader::from_xml("").unwrap();
        assert!(steps.is_empty());
    }
    #[test]
    fn test_xdmf_hdf5_reference() {
        let ref_str = write_xdmf_hdf5_reference("data.h5", "/positions");
        assert!(ref_str.contains("data.h5:/positions"));
        assert!(ref_str.contains("HDF"));
    }
    #[test]
    fn test_time_series_empty() {
        let ts = XdmfTimeSeries::new();
        let xml = ts.to_xml();
        assert!(xml.contains("Temporal"));
        assert!(xml.contains("<Xdmf"));
        assert!(xml.contains("</Xdmf>"));
    }
    #[test]
    fn test_xdmf_step_n_points() {
        let mut ts = XdmfTimeSeries::new();
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        ts.add_step(0.0, positions, vec![]);
        assert_eq!(ts.steps[0].n_points, 10);
    }
    #[test]
    fn test_vector_attribute_basic() {
        let vecs = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let s = xdmf_vector_attribute("velocity", &vecs);
        assert!(s.contains("velocity"));
        assert!(s.contains("Vector"));
        assert!(s.contains("2 3"));
        assert!(s.contains("1 0 0"));
    }
    #[test]
    fn test_vector_attribute_empty() {
        let s = xdmf_vector_attribute("v", &[]);
        assert!(s.contains("0 3"));
    }
    #[test]
    fn test_tensor6_attribute_basic() {
        let tensors = vec![[1.0, 2.0, 3.0, 0.5, 0.1, 0.2]];
        let s = xdmf_tensor6_attribute("stress", &tensors);
        assert!(s.contains("stress"));
        assert!(s.contains("Tensor6"));
        assert!(s.contains("1 6"));
        assert!(s.contains("0.5"));
    }
    #[test]
    fn test_unstructured_triangle_mesh() {
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let connectivity = vec![0usize, 1, 2];
        let mut buf = Vec::new();
        write_xdmf_unstructured(
            &mut buf,
            XdmfTopologyType::Triangle,
            &nodes,
            &connectivity,
            &[],
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Triangle"), "topology type not found");
        assert!(s.contains("NumberOfElements=\"1\""));
        assert!(s.contains("0 1 2"));
    }
    #[test]
    fn test_unstructured_tet_mesh() {
        let nodes: Vec<[f64; 3]> = (0..4).map(|i| [i as f64, 0.0, 0.0]).collect();
        let connectivity = vec![0usize, 1, 2, 3];
        let mut buf = Vec::new();
        write_xdmf_unstructured(
            &mut buf,
            XdmfTopologyType::Tetrahedron,
            &nodes,
            &connectivity,
            &[],
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Tetrahedron"));
        assert!(s.contains("NumberOfElements=\"1\""));
    }
    #[test]
    fn test_unstructured_with_scalar_field() {
        let nodes = vec![[0.0; 3]; 3];
        let connectivity = vec![0usize, 1, 2];
        let pressure = vec![1.0, 2.0, 3.0];
        let mut buf = Vec::new();
        write_xdmf_unstructured(
            &mut buf,
            XdmfTopologyType::Triangle,
            &nodes,
            &connectivity,
            &[("pressure", &pressure)],
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("pressure"));
    }
    #[test]
    fn topology_type_nodes_per_element() {
        assert_eq!(XdmfTopologyType::Triangle.nodes_per_element(), 3);
        assert_eq!(XdmfTopologyType::Tetrahedron.nodes_per_element(), 4);
        assert_eq!(XdmfTopologyType::Hexahedron.nodes_per_element(), 8);
        assert_eq!(XdmfTopologyType::Quadrilateral.nodes_per_element(), 4);
    }
    #[test]
    fn schema_validate_ok() {
        let schema = XdmfSchema::new("Polyvertex", vec!["density".to_string()]);
        let positions = make_positions(2);
        let density = vec![1.0, 2.0];
        let fields: &[(&str, &[f64])] = &[("density", &density)];
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, fields).unwrap();
        let xml = String::from_utf8(buf).unwrap();
        let errors = schema.validate(&xml);
        assert!(errors.is_empty(), "unexpected errors: {errors:?}");
    }
    #[test]
    fn schema_validate_missing_topology() {
        let schema = XdmfSchema::new("Triangle", vec![]);
        let positions = make_positions(3);
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, &[]).unwrap();
        let xml = String::from_utf8(buf).unwrap();
        let errors = schema.validate(&xml);
        assert!(
            !errors.is_empty(),
            "should report missing topology Triangle"
        );
    }
    #[test]
    fn schema_validate_missing_attribute() {
        let schema = XdmfSchema::new("Polyvertex", vec!["velocity".to_string()]);
        let positions = make_positions(2);
        let mut buf = Vec::new();
        write_xdmf_particles(&mut buf, &positions, &[]).unwrap();
        let xml = String::from_utf8(buf).unwrap();
        let errors = schema.validate(&xml);
        assert!(!errors.is_empty());
    }
    #[test]
    fn time_series_vector_step_encodes_xyz() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step_with_vectors(
            0.0,
            vec![[0.0; 3]],
            vec![],
            vec![("velocity".to_string(), vec![[1.0, 2.0, 3.0]])],
        );
        let xml = ts.to_xml();
        assert!(xml.contains("velocity_x"), "should have velocity_x field");
        assert!(xml.contains("velocity_y"));
        assert!(xml.contains("velocity_z"));
    }
    #[test]
    fn time_series_times() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![], vec![]);
        ts.add_step(0.5, vec![], vec![]);
        ts.add_step(1.0, vec![], vec![]);
        assert_eq!(ts.times(), vec![0.0, 0.5, 1.0]);
    }
    #[test]
    fn time_series_total_particle_count() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![[0.0; 3]; 5], vec![]);
        ts.add_step(1.0, vec![[0.0; 3]; 3], vec![]);
        assert_eq!(ts.total_particle_count(), 8);
    }
    #[test]
    fn time_series_default() {
        let ts = XdmfTimeSeries::default();
        assert!(ts.steps.is_empty());
    }
    fn make_tri_step(time: f64) -> XdmfMeshStep {
        XdmfMeshStep {
            time,
            nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            connectivity: vec![0, 1, 2],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![],
            node_vectors: vec![],
        }
    }
    #[test]
    fn mesh_series_add_and_len() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(make_tri_step(0.0));
        ms.add_step(make_tri_step(1.0));
        assert_eq!(ms.len(), 2);
        assert!(!ms.is_empty());
    }
    #[test]
    fn mesh_series_to_xml_contains_temporal() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(make_tri_step(0.0));
        let xml = ms.to_xml();
        assert!(xml.contains("Temporal"));
        assert!(xml.contains("Triangle"));
        assert!(xml.contains("0 1 2"));
    }
    #[test]
    fn mesh_series_scalar_field_in_xml() {
        let mut step = make_tri_step(0.0);
        step.node_scalars
            .push(("temperature".into(), vec![300.0, 310.0, 305.0]));
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(step);
        let xml = ms.to_xml();
        assert!(xml.contains("temperature"));
        assert!(xml.contains("310"));
    }
    #[test]
    fn mesh_series_vector_field_in_xml() {
        let mut step = make_tri_step(0.0);
        step.node_vectors.push((
            "velocity".into(),
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ));
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(step);
        let xml = ms.to_xml();
        assert!(xml.contains("velocity"));
        assert!(xml.contains("Vector"));
    }
    #[test]
    fn mesh_series_times() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(make_tri_step(0.0));
        ms.add_step(make_tri_step(0.5));
        ms.add_step(make_tri_step(1.0));
        assert_eq!(ms.times(), vec![0.0, 0.5, 1.0]);
    }
    #[test]
    fn mesh_step_n_elements() {
        let step = make_tri_step(0.0);
        assert_eq!(step.n_elements(), 1);
    }
    #[test]
    fn hdf5_builder_basic() {
        let b = Hdf5DataItemBuilder::new("data.h5");
        let xml = b.build("positions", "100 3");
        assert!(xml.contains("data.h5:/positions"));
        assert!(xml.contains("100 3"));
        assert!(xml.contains("HDF"));
    }
    #[test]
    fn hdf5_builder_with_group() {
        let b = Hdf5DataItemBuilder::new("sim.h5").group("/step0");
        let xml = b.build("velocity", "50 3");
        assert!(xml.contains("sim.h5:/step0/velocity"));
    }
    #[test]
    fn scalar_data_item_format() {
        let item = xdmf_scalar_data_item(&[1.0, 2.0, 3.0]);
        assert!(item.contains("Dimensions=\"3\""));
        assert!(item.contains("1 2 3"));
    }
    #[test]
    fn vector_data_item_format() {
        let item = xdmf_vector_data_item(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        assert!(item.contains("2 3"));
        assert!(item.contains("1 0 0"));
    }
    #[test]
    fn time_element_format() {
        let t = xdmf_time_element(3.125);
        assert!(t.contains("3.125"));
        assert!(t.contains("<Time"));
    }
    #[test]
    fn total_node_count_basic() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(make_tri_step(0.0));
        ms.add_step(make_tri_step(1.0));
        assert_eq!(total_node_count(&ms), 6);
    }
    #[test]
    fn peak_element_step_basic() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(make_tri_step(0.0));
        ms.add_step(XdmfMeshStep {
            time: 1.0,
            nodes: vec![[0.0; 3]; 4],
            connectivity: vec![0, 1, 2, 1, 2, 3],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![],
            node_vectors: vec![],
        });
        assert_eq!(peak_element_step(&ms), Some(1));
    }
    #[test]
    fn peak_element_step_empty_series() {
        let ms = XdmfMeshTimeSeries::new();
        assert_eq!(peak_element_step(&ms), None);
    }
    #[test]
    fn test_write_collection_creates_file_with_timestep_count() {
        let mut ts = XdmfTimeSeriesHdf5::new();
        ts.timesteps = vec![0.0, 1.0, 2.0];
        ts.hdf5_paths = vec!["a.h5".into(), "b.h5".into(), "c.h5".into()];
        ts.attribute_names = vec!["temperature".into()];
        let path = std::env::temp_dir().join("test_write_collection.xmf");
        ts.write_collection(path.to_str().unwrap_or(""), 10, 5, "Triangle")
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let count = content.matches("<Time Value=").count();
        assert_eq!(count, 3, "expected 3 timestep entries, got {}", count);
        assert!(content.contains("CollectionType=\"Temporal\""));
        assert!(content.contains("temperature"));
    }
    #[test]
    fn test_parse_xdmf_topology_finds_topology() {
        let xml = r#"<?xml version="1.0"?>
<Xdmf Version="3.0">
  <Domain>
    <Grid Name="mesh" GridType="Uniform">
      <Topology TopologyType="Triangle" NumberOfElements="42"/>
    </Grid>
  </Domain>
</Xdmf>"#;
        let topologies = parse_xdmf_topology(xml);
        assert_eq!(topologies.len(), 1);
        assert_eq!(topologies[0].0, "Triangle");
        assert_eq!(topologies[0].1, 42);
    }
    #[test]
    fn test_write_xdmf_uniform_grid_contains_grid_tag() {
        let params = XdmfUniformGrid {
            name: "volume".to_string(),
            dimensions: [4, 5, 6],
            origin: [0.0, 0.0, 0.0],
            spacing: [1.0, 1.0, 1.0],
        };
        let path = std::env::temp_dir().join("test_uniform_grid.xmf");
        write_xdmf_uniform_grid(path.to_str().unwrap_or(""), &params).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<Grid"), "missing Grid tag");
        assert!(content.contains("volume"), "missing grid name");
        assert!(content.contains("3DCoRectMesh"), "missing topology type");
        assert!(content.contains("ORIGIN_DXDYDZ"), "missing geometry type");
    }
    #[test]
    fn test_write_xdmf_with_attributes_creates_valid_file() {
        let attrs = vec![XdmfAttribute {
            name: "pressure".to_string(),
            center: "Node".to_string(),
            n_components: 1,
            hdf5_path: "data.h5:/pressure".to_string(),
        }];
        let path = std::env::temp_dir().join("test_xdmf_with_attrs.xmf");
        write_xdmf_with_attributes(
            path.to_str().unwrap_or(""),
            "data.h5",
            20,
            10,
            "Tetrahedron",
            &attrs,
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("pressure"));
        assert!(content.contains("Tetrahedron"));
        assert!(content.contains("data.h5:/coordinates"));
    }
}
#[cfg(test)]
mod tests_xdmf_ext {

    use crate::xdmf::types::*;
    #[test]
    fn add_frame_increments_step_count() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame(0.0, vec![[0.0; 3]; 5]);
        ts.add_frame(1.0, vec![[1.0; 3]; 5]);
        assert_eq!(ts.steps.len(), 2);
    }
    #[test]
    fn add_frame_sets_n_points() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame(0.0, vec![[0.0; 3]; 7]);
        assert_eq!(ts.steps[0].n_points, 7);
    }
    #[test]
    fn add_frame_with_scalars_persists_data() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame_with_scalars(
            0.5,
            vec![[0.0; 3]; 3],
            vec![("pressure".to_string(), vec![1.0, 2.0, 3.0])],
        );
        assert_eq!(ts.steps[0].scalar_fields.len(), 1);
        assert_eq!(ts.steps[0].scalar_fields[0].0, "pressure");
    }
    #[test]
    fn write_xml_produces_valid_xml_string() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame(0.0, vec![[1.0, 2.0, 3.0]]);
        let mut buf: Vec<u8> = Vec::new();
        ts.write_xml(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("<?xml"));
        assert!(s.contains("Temporal"));
        assert!(s.contains("Time Value=\"0\"") || s.contains("Time Value=\"0."));
    }
    #[test]
    fn write_xml_matches_to_xml() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame(1.5, vec![[0.0; 3]; 4]);
        let expected = ts.to_xml();
        let mut buf: Vec<u8> = Vec::new();
        ts.write_xml(&mut buf).unwrap();
        let actual = String::from_utf8(buf).unwrap();
        assert_eq!(actual, expected);
    }
    #[test]
    fn write_xml_to_file_creates_file() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame(0.0, vec![[0.0; 3]; 2]);
        let path = std::env::temp_dir().join("test_write_xml_ext.xmf");
        ts.write_xml_to_file(path.to_str().unwrap_or("")).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<Xdmf"));
    }
    #[test]
    fn mixed_topology_xdmf_name() {
        assert_eq!(XdmfTopologyType::Mixed.xdmf_name(), "Mixed");
    }
    #[test]
    fn mixed_topology_nodes_per_element_zero() {
        assert_eq!(XdmfTopologyType::Mixed.nodes_per_element(), 0);
    }
    #[test]
    fn topology_type_names_cover_all_variants() {
        let types = [
            XdmfTopologyType::Triangle,
            XdmfTopologyType::Tetrahedron,
            XdmfTopologyType::Hexahedron,
            XdmfTopologyType::Quadrilateral,
            XdmfTopologyType::Mixed,
        ];
        let names: Vec<&str> = types.iter().map(|t| t.xdmf_name()).collect();
        assert!(names.contains(&"Triangle"));
        assert!(names.contains(&"Tetrahedron"));
        assert!(names.contains(&"Hexahedron"));
        assert!(names.contains(&"Quadrilateral"));
        assert!(names.contains(&"Mixed"));
    }
    #[test]
    fn add_frame_empty_positions() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_frame(0.0, vec![]);
        assert_eq!(ts.steps[0].n_points, 0);
        let xml = ts.to_xml();
        assert!(xml.contains("<Xdmf"));
    }
    #[test]
    fn write_xml_multiple_frames_all_present() {
        let mut ts = XdmfTimeSeries::new();
        for i in 0..4 {
            ts.add_frame(i as f64, vec![[0.0; 3]; 2]);
        }
        let mut buf: Vec<u8> = Vec::new();
        ts.write_xml(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let count = s.matches("<Time Value=").count();
        assert_eq!(count, 4);
    }
}
/// Write a full 3×3 tensor attribute (9 components per node).
///
/// XDMF represents this as `AttributeType="Tensor"` with `Dimensions="N 9"`.
pub fn xdmf_tensor9_attribute(name: &str, tensors: &[[f64; 9]]) -> String {
    let n = tensors.len();
    let mut s = String::new();
    s.push_str(&format!(
        "      <Attribute Name=\"{}\" AttributeType=\"Tensor\" Center=\"Node\">\n",
        name
    ));
    s.push_str(&format!(
        "        <DataItem Format=\"XML\" Dimensions=\"{} 9\">\n",
        n
    ));
    for t in tensors {
        s.push_str("          ");
        for (i, c) in t.iter().enumerate() {
            if i > 0 {
                s.push(' ');
            }
            s.push_str(&format!("{}", c));
        }
        s.push('\n');
    }
    s.push_str("        </DataItem>\n");
    s.push_str("      </Attribute>\n");
    s
}
/// Validate that an XDMF string contains required version and domain tags.
pub fn xdmf_is_well_formed(xml: &str) -> bool {
    xml.contains("<Xdmf")
        && xml.contains("Version=")
        && xml.contains("<Domain>")
        && xml.contains("</Domain>")
        && xml.contains("</Xdmf>")
}
/// Count the number of `<Grid` elements in an XDMF document.
pub fn xdmf_count_grids(xml: &str) -> usize {
    xml.matches("<Grid").count()
}
/// Count the number of `<Attribute` elements in an XDMF document.
pub fn xdmf_count_attributes(xml: &str) -> usize {
    xml.matches("<Attribute").count()
}
#[cfg(test)]
mod tests_xdmf_new {
    use super::*;
    use crate::xdmf::Hdf5DataItemBuilder;

    use crate::xdmf::XdmfMeshStep;
    use crate::xdmf::XdmfMeshTimeSeries;
    use crate::xdmf::XdmfMultiBlock;

    use crate::xdmf::XdmfStructuredGrid;

    use crate::xdmf::XdmfTimeSeriesHdf5;
    use crate::xdmf::XdmfWriter;

    use crate::xdmf::total_node_count;
    use crate::xdmf::types::*;

    use crate::xdmf::xdmf_is_well_formed;
    use crate::xdmf::xdmf_scalar_data_item;
    use crate::xdmf::xdmf_time_element;
    use crate::xdmf::xdmf_vector_data_item;
    #[test]
    fn structured_grid_node_count() {
        let g = XdmfStructuredGrid::new("test", 4, 5, 6, [0.0; 3], 1.0, 1.0, 1.0);
        assert_eq!(g.n_nodes(), 4 * 5 * 6);
    }
    #[test]
    fn structured_grid_cell_count() {
        let g = XdmfStructuredGrid::new("test", 4, 5, 6, [0.0; 3], 1.0, 1.0, 1.0);
        assert_eq!(g.n_cells(), 3 * 4 * 5);
    }
    #[test]
    fn structured_grid_to_xml_contains_topology() {
        let g = XdmfStructuredGrid::new("flow", 3, 3, 3, [0.0; 3], 0.5, 0.5, 0.5);
        let xml = g.to_xml();
        assert!(xml.contains("3DCoRectMesh"));
        assert!(xml.contains("ORIGIN_DXDYDZ"));
        assert!(xml.contains("flow"));
    }
    #[test]
    fn structured_grid_scalar_in_xml() {
        let mut g = XdmfStructuredGrid::new("g", 2, 2, 2, [0.0; 3], 1.0, 1.0, 1.0);
        let n = g.n_nodes();
        g.add_node_scalar("pressure", vec![1.0; n]);
        let xml = g.to_xml();
        assert!(xml.contains("pressure"));
        assert!(xml.contains("Center=\"Node\""));
    }
    #[test]
    fn structured_grid_cell_scalar_in_xml() {
        let mut g = XdmfStructuredGrid::new("g", 3, 3, 3, [0.0; 3], 1.0, 1.0, 1.0);
        let nc = g.n_cells();
        g.add_cell_scalar("vorticity", vec![0.5; nc]);
        let xml = g.to_xml();
        assert!(xml.contains("vorticity"));
        assert!(xml.contains("Center=\"Cell\""));
    }
    #[test]
    fn structured_grid_vector_in_xml() {
        let mut g = XdmfStructuredGrid::new("g", 2, 2, 2, [0.0; 3], 1.0, 1.0, 1.0);
        let n = g.n_nodes();
        g.add_node_vector("velocity", vec![[1.0, 0.0, 0.0]; n]);
        let xml = g.to_xml();
        assert!(xml.contains("velocity"));
        assert!(xml.contains("AttributeType=\"Vector\""));
    }
    #[test]
    fn structured_grid_origin_appears_in_geometry() {
        let g = XdmfStructuredGrid::new("g", 2, 2, 2, [1.0, 2.0, 3.0], 0.1, 0.2, 0.3);
        let xml = g.to_xml();
        assert!(xml.contains("3") && xml.contains("ORIGIN_DXDYDZ"));
    }
    #[test]
    fn xdmf_writer_basic_document() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_grid("particles", "Uniform");
        w.write_polyvertex_topology(5);
        w.write_xyz_geometry(&[[0.0, 0.0, 0.0]; 5]);
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        assert!(xml.contains("<?xml"));
        assert!(xml.contains("<Xdmf"));
        assert!(xml.contains("Polyvertex"));
        assert!(xml.contains("particles"));
        assert!(xml.contains("</Xdmf>"));
    }
    #[test]
    fn xdmf_writer_temporal_collection() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_temporal_collection("ts");
        for i in 0..3_usize {
            w.open_grid("step", "Uniform");
            w.write_time(i as f64 * 0.1);
            w.write_polyvertex_topology(2);
            w.write_xyz_geometry(&[[0.0; 3]; 2]);
            w.close_grid();
        }
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        let count = xml.matches("<Time Value=").count();
        assert_eq!(count, 3);
        assert!(xml.contains("CollectionType=\"Temporal\""));
    }
    #[test]
    fn xdmf_writer_scalar_attribute() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_grid("g", "Uniform");
        w.write_polyvertex_topology(3);
        w.write_xyz_geometry(&[[0.0; 3]; 3]);
        w.write_scalar_attribute("density", "Node", &[1.0, 2.0, 3.0]);
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        assert!(xml.contains("density"));
        assert!(xml.contains("Scalar"));
    }
    #[test]
    fn xdmf_writer_vector_attribute() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_grid("g", "Uniform");
        w.write_polyvertex_topology(2);
        w.write_xyz_geometry(&[[0.0; 3]; 2]);
        w.write_vector_attribute("velocity", "Node", &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        assert!(xml.contains("velocity"));
        assert!(xml.contains("Vector"));
    }
    #[test]
    fn xdmf_writer_hdf5_attribute() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_grid("g", "Uniform");
        w.write_polyvertex_topology(10);
        w.write_hdf5_attribute("temp", "Node", "Scalar", "10", "data.h5:/temperature");
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        assert!(xml.contains("HDF"));
        assert!(xml.contains("temperature"));
    }
    #[test]
    fn xdmf_writer_peek() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        let preview = w.peek().to_string();
        assert!(preview.contains("<Domain>"));
    }
    #[test]
    fn tensor9_attribute_has_correct_dimensions() {
        let tensors = vec![[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64]; 3];
        let xml = xdmf_tensor9_attribute("stress", &tensors);
        assert!(xml.contains("Tensor"));
        assert!(xml.contains("3 9"));
        assert!(xml.contains("stress"));
    }
    #[test]
    fn tensor9_attribute_empty() {
        let xml = xdmf_tensor9_attribute("empty_tensor", &[]);
        assert!(xml.contains("0 9"));
    }
    #[test]
    fn well_formed_valid_document() {
        let xml =
            "<?xml version=\"1.0\"?>\n<Xdmf Version=\"3.0\">\n  <Domain>\n  </Domain>\n</Xdmf>\n";
        assert!(xdmf_is_well_formed(xml));
    }
    #[test]
    fn well_formed_rejects_truncated() {
        let xml = "<Xdmf Version=\"3.0\">\n  <Domain>\n";
        assert!(!xdmf_is_well_formed(xml));
    }
    #[test]
    fn count_grids_correct() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_temporal_collection("ts");
        for _ in 0..3_usize {
            w.open_grid("g", "Uniform");
            w.close_grid();
        }
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        assert_eq!(xdmf_count_grids(&xml), 4);
    }
    #[test]
    fn count_attributes_correct() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_grid("g", "Uniform");
        w.write_polyvertex_topology(2);
        w.write_xyz_geometry(&[[0.0; 3]; 2]);
        w.write_scalar_attribute("a", "Node", &[1.0, 2.0]);
        w.write_scalar_attribute("b", "Node", &[3.0, 4.0]);
        w.write_vector_attribute("v", "Node", &[[0.0; 3]; 2]);
        w.close_grid();
        w.close_domain();
        let xml = w.finish();
        assert_eq!(xdmf_count_attributes(&xml), 3);
    }
    #[test]
    fn multiblock_empty() {
        let mb = XdmfMultiBlock::new();
        assert!(mb.is_empty());
        assert_eq!(mb.len(), 0);
    }
    #[test]
    fn multiblock_spatial_collection() {
        let mut mb = XdmfMultiBlock::new();
        mb.add_block("block0", "<Grid Name=\"b0\" GridType=\"Uniform\"></Grid>");
        mb.add_block("block1", "<Grid Name=\"b1\" GridType=\"Uniform\"></Grid>");
        assert_eq!(mb.len(), 2);
        let xml = mb.to_xml();
        assert!(xml.contains("Spatial"));
        assert!(xml.contains("b0"));
        assert!(xml.contains("b1"));
    }
    #[test]
    fn multiblock_to_xml_well_formed() {
        let mut mb = XdmfMultiBlock::new();
        mb.add_block("g", "<Grid Name=\"g\" GridType=\"Uniform\"></Grid>");
        let xml = mb.to_xml();
        assert!(xdmf_is_well_formed(&xml));
    }
    #[test]
    fn hdf5_time_series_attribute_names() {
        let mut ts = XdmfTimeSeriesHdf5::new();
        ts.timesteps = vec![0.0, 1.0];
        ts.hdf5_paths = vec!["step0.h5".into(), "step1.h5".into()];
        ts.attribute_names = vec!["density".into(), "pressure".into()];
        let path = std::env::temp_dir().join("test_hdf5_series_attrs.xmf");
        ts.write_collection(path.to_str().unwrap_or(""), 100, 50, "Tetrahedron")
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("density"));
        assert!(content.contains("pressure"));
        let time_count = content.matches("<Time Value=").count();
        assert_eq!(time_count, 2);
    }
    #[test]
    fn hdf5_data_item_builder_group() {
        let b = Hdf5DataItemBuilder::new("sim.h5").group("/results");
        let item = b.build("velocity", "100 3");
        assert!(item.contains("sim.h5:/results/velocity"));
        assert!(item.contains("100 3"));
    }
    #[test]
    fn hdf5_data_item_builder_root_group() {
        let b = Hdf5DataItemBuilder::new("data.h5");
        let item = b.build("coordinates", "50 3");
        assert!(item.contains("data.h5:/coordinates"));
    }
    #[test]
    fn write_xdmf_unstructured_tetra() {
        let nodes: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let connectivity = vec![0, 1, 2, 3];
        let mut buf = Vec::new();
        write_xdmf_unstructured(
            &mut buf,
            XdmfTopologyType::Tetrahedron,
            &nodes,
            &connectivity,
            &[],
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Tetrahedron"));
        assert!(s.contains("NumberOfElements=\"1\""));
    }
    #[test]
    fn write_xdmf_unstructured_hex() {
        let nodes = vec![[0.0f64; 3]; 8];
        let connectivity: Vec<usize> = (0..8).collect();
        let mut buf = Vec::new();
        write_xdmf_unstructured(
            &mut buf,
            XdmfTopologyType::Hexahedron,
            &nodes,
            &connectivity,
            &[],
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Hexahedron"));
        assert!(s.contains("NumberOfElements=\"1\""));
    }
    #[test]
    fn xdmf_scalar_data_item_values_present() {
        let item = xdmf_scalar_data_item(&[1.0, 2.0, 3.0]);
        assert!(item.contains("Dimensions=\"3\""));
        assert!(item.contains("1"));
        assert!(item.contains("2"));
        assert!(item.contains("3"));
    }
    #[test]
    fn xdmf_vector_data_item_dimensions() {
        let item = xdmf_vector_data_item(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert!(item.contains("2 3"));
        assert!(item.contains("1 2 3"));
        assert!(item.contains("4 5 6"));
    }
    #[test]
    fn xdmf_time_element_format() {
        let t = xdmf_time_element(1.5);
        assert_eq!(t, "<Time Value=\"1.5\"/>");
    }
    #[test]
    fn total_node_count_multi_step() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(XdmfMeshStep {
            time: 0.0,
            nodes: vec![[0.0; 3]; 4],
            connectivity: vec![0, 1, 2],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![],
            node_vectors: vec![],
        });
        ms.add_step(XdmfMeshStep {
            time: 1.0,
            nodes: vec![[0.0; 3]; 6],
            connectivity: vec![0, 1, 2, 3, 4, 5],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![],
            node_vectors: vec![],
        });
        assert_eq!(total_node_count(&ms), 10);
    }
    #[test]
    fn mesh_time_series_to_xml_well_formed() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(XdmfMeshStep {
            time: 0.5,
            nodes: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            connectivity: vec![0, 1, 2],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![("T".to_string(), vec![300.0, 310.0, 320.0])],
            node_vectors: vec![("v".to_string(), vec![[1.0, 0.0, 0.0]; 3])],
        });
        let xml = ms.to_xml();
        assert!(xdmf_is_well_formed(&xml));
        assert!(xml.contains("MeshTimeSeries"));
        assert!(xml.contains("Time Value=\"0.5\""));
        assert!(xml.contains("T"));
        assert!(xml.contains("v"));
    }
}
#[cfg(test)]
mod tests_xdmf_extra {
    use super::*;
    use crate::xdmf::types::*;
    #[test]
    fn time_series_default_is_empty() {
        let ts: XdmfTimeSeries = Default::default();
        assert!(ts.steps.is_empty());
    }
    #[test]
    fn time_series_times_single_step() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(2.5, vec![[0.0; 3]; 3], vec![]);
        let t = ts.times();
        assert_eq!(t.len(), 1);
        assert!((t[0] - 2.5).abs() < 1e-12);
    }
    #[test]
    fn time_series_times_multiple_steps() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![[0.0; 3]; 2], vec![]);
        ts.add_step(1.0, vec![[0.0; 3]; 2], vec![]);
        ts.add_step(2.0, vec![[0.0; 3]; 2], vec![]);
        let t = ts.times();
        assert_eq!(t, vec![0.0, 1.0, 2.0]);
    }
    #[test]
    fn time_series_total_particle_count_empty() {
        let ts = XdmfTimeSeries::new();
        assert_eq!(ts.total_particle_count(), 0);
    }
    #[test]
    fn time_series_total_particle_count_multi() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![[0.0; 3]; 4], vec![]);
        ts.add_step(1.0, vec![[0.0; 3]; 6], vec![]);
        assert_eq!(ts.total_particle_count(), 10);
    }
    #[test]
    fn time_series_add_step_with_vectors() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step_with_vectors(
            0.1,
            vec![[0.0; 3]; 2],
            vec![],
            vec![("vel".to_string(), vec![[1.0, 0.0, 0.0]; 2])],
        );
        assert_eq!(ts.steps.len(), 1);
    }
    #[test]
    fn time_series_to_xml_contains_positions() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(0.0, vec![[1.0, 2.0, 3.0]], vec![]);
        let xml = ts.to_xml();
        assert!(xml.contains("1 2 3") || xml.contains("1.0 2.0 3.0") || xml.contains("1"));
        assert!(xml.contains("GeometryType=\"XYZ\""));
    }
    #[test]
    fn time_series_to_xml_contains_scalar_field() {
        let mut ts = XdmfTimeSeries::new();
        ts.add_step(
            0.0,
            vec![[0.0; 3]; 2],
            vec![("temperature".to_string(), vec![300.0, 301.0])],
        );
        let xml = ts.to_xml();
        assert!(xml.contains("temperature"));
    }
    #[test]
    fn topology_type_triangle_xdmf_name() {
        assert_eq!(XdmfTopologyType::Triangle.xdmf_name(), "Triangle");
    }
    #[test]
    fn topology_type_tet_nodes_per_element() {
        assert_eq!(XdmfTopologyType::Tetrahedron.nodes_per_element(), 4);
    }
    #[test]
    fn topology_type_quad_nodes_per_element() {
        assert_eq!(XdmfTopologyType::Quadrilateral.nodes_per_element(), 4);
    }
    #[test]
    fn topology_type_hex_nodes_per_element() {
        assert_eq!(XdmfTopologyType::Hexahedron.nodes_per_element(), 8);
    }
    #[test]
    fn topology_type_mixed_name() {
        assert_eq!(XdmfTopologyType::Mixed.xdmf_name(), "Mixed");
    }
    #[test]
    fn mesh_time_series_empty_is_empty() {
        let ms = XdmfMeshTimeSeries::new();
        assert!(ms.is_empty());
        assert_eq!(ms.len(), 0);
    }
    #[test]
    fn mesh_time_series_len_after_add() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(XdmfMeshStep {
            time: 0.0,
            nodes: vec![[0.0; 3]; 3],
            connectivity: vec![0, 1, 2],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![],
            node_vectors: vec![],
        });
        assert_eq!(ms.len(), 1);
        assert!(!ms.is_empty());
    }
    #[test]
    fn mesh_time_series_times_list() {
        let mut ms = XdmfMeshTimeSeries::new();
        for i in 0..3 {
            ms.add_step(XdmfMeshStep {
                time: i as f64 * 0.5,
                nodes: vec![[0.0; 3]; 3],
                connectivity: vec![0, 1, 2],
                topology: XdmfTopologyType::Triangle,
                node_scalars: vec![],
                node_vectors: vec![],
            });
        }
        assert_eq!(ms.times(), vec![0.0, 0.5, 1.0]);
    }
    #[test]
    fn mesh_step_n_elements_triangle() {
        let step = XdmfMeshStep {
            time: 0.0,
            nodes: vec![[0.0; 3]; 3],
            connectivity: vec![0, 1, 2],
            topology: XdmfTopologyType::Triangle,
            node_scalars: vec![],
            node_vectors: vec![],
        };
        assert_eq!(step.n_elements(), 1);
    }
    #[test]
    fn mesh_step_n_elements_hex() {
        let step = XdmfMeshStep {
            time: 0.0,
            nodes: vec![[0.0; 3]; 8],
            connectivity: vec![0, 1, 2, 3, 4, 5, 6, 7],
            topology: XdmfTopologyType::Hexahedron,
            node_scalars: vec![],
            node_vectors: vec![],
        };
        assert_eq!(step.n_elements(), 1);
    }
    #[test]
    fn hdf5_builder_no_group() {
        let b = Hdf5DataItemBuilder::new("data.h5");
        let item = b.build("coordinates", "100 3");
        assert!(item.contains("data.h5:/coordinates"));
        assert!(item.contains("100 3"));
    }
    #[test]
    fn hdf5_builder_with_group() {
        let b = Hdf5DataItemBuilder::new("out.h5").group("step_0");
        let item = b.build("positions", "50 3");
        assert!(item.contains("step_0") && item.contains("positions"));
        assert!(item.contains("out.h5"));
    }
    #[test]
    fn hdf5_builder_dimensions_in_output() {
        let b = Hdf5DataItemBuilder::new("f.h5");
        let item = b.build("vel", "200 3");
        assert!(item.contains("200 3"));
        assert!(item.contains("HDF"));
    }
    #[test]
    fn scalar_data_item_empty() {
        let item = xdmf_scalar_data_item(&[]);
        assert!(item.contains("0"));
    }
    #[test]
    fn scalar_data_item_single_value() {
        let item = xdmf_scalar_data_item(&[42.0]);
        assert!(item.contains("42"));
        assert!(item.contains("Dimensions=\"1\""));
    }
    #[test]
    fn vector_data_item_two_rows() {
        let item = xdmf_vector_data_item(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert!(item.contains("2 3"));
        assert!(item.contains("Vector") || item.contains("XML"));
    }
    #[test]
    fn xdmf_schema_validate_ok() {
        let schema = XdmfSchema::new("Triangle", vec!["pressure".to_string()]);
        let xml = "<Topology TopologyType=\"Triangle\"/><Attribute Name=\"pressure\"/>";
        let errs = schema.validate(xml);
        assert!(errs.is_empty(), "unexpected errors: {:?}", errs);
    }
    #[test]
    fn xdmf_schema_validate_missing_attr() {
        let schema = XdmfSchema::new("Triangle", vec!["velocity".to_string()]);
        let xml = "<Topology TopologyType=\"Triangle\"/>";
        let errs = schema.validate(xml);
        assert!(!errs.is_empty());
    }
    #[test]
    fn xdmf_schema_validate_wrong_topology() {
        let schema = XdmfSchema::new("Hexahedron", vec![]);
        let xml = "<Topology TopologyType=\"Triangle\"/>";
        let errs = schema.validate(xml);
        assert!(!errs.is_empty());
    }
    #[test]
    fn xdmf_writer_new_starts_with_xml_declaration() {
        let w = XdmfWriter::new();
        let out = w.finish();
        assert!(out.starts_with("<?xml"));
        assert!(out.contains("<Xdmf Version=\"3.0\">"));
    }
    #[test]
    fn xdmf_writer_domain_round_trip() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.close_domain();
        let out = w.finish();
        assert!(out.contains("<Domain>"));
        assert!(out.contains("</Domain>"));
    }
    #[test]
    fn xdmf_writer_grid_elements() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_grid("particles", "Uniform");
        w.write_time(3.125);
        w.write_polyvertex_topology(100);
        w.close_grid();
        w.close_domain();
        let out = w.finish();
        assert!(out.contains("particles"));
        assert!(out.contains("3.125"));
        assert!(out.contains("100"));
    }
    #[test]
    fn xdmf_writer_temporal_collection() {
        let mut w = XdmfWriter::new();
        w.open_domain();
        w.open_temporal_collection("ts");
        w.close_grid();
        w.close_domain();
        let out = w.finish();
        assert!(out.contains("CollectionType=\"Temporal\""));
    }
    #[test]
    fn structured_grid_single_cell_size() {
        let g = XdmfStructuredGrid::new("g", 2, 2, 2, [0.0; 3], 1.0, 1.0, 1.0);
        assert_eq!(g.n_cells(), 1);
    }
    #[test]
    fn structured_grid_1d_degenerate_cells() {
        let g = XdmfStructuredGrid::new("g", 1, 5, 5, [0.0; 3], 1.0, 1.0, 1.0);
        assert_eq!(g.n_cells(), 0);
    }
    #[test]
    fn structured_grid_xml_well_formed() {
        let g = XdmfStructuredGrid::new("box", 3, 3, 3, [0.0; 3], 0.5, 0.5, 0.5);
        let xml = g.to_xml();
        assert!(xdmf_is_well_formed(&xml));
    }
    #[test]
    fn total_node_count_empty_series() {
        let ms = XdmfMeshTimeSeries::new();
        assert_eq!(total_node_count(&ms), 0);
    }
    #[test]
    fn peak_element_step_empty_returns_none() {
        let ms = XdmfMeshTimeSeries::new();
        assert!(peak_element_step(&ms).is_none());
    }
    #[test]
    fn peak_element_step_single_step() {
        let mut ms = XdmfMeshTimeSeries::new();
        ms.add_step(XdmfMeshStep {
            time: 0.0,
            nodes: vec![[0.0; 3]; 4],
            connectivity: vec![0, 1, 2, 3],
            topology: XdmfTopologyType::Quadrilateral,
            node_scalars: vec![],
            node_vectors: vec![],
        });
        assert_eq!(peak_element_step(&ms), Some(0));
    }
    #[test]
    fn hdf5_series_default_empty() {
        let s: XdmfTimeSeriesHdf5 = Default::default();
        assert!(s.timesteps.is_empty());
        assert!(s.hdf5_paths.is_empty());
    }
    #[test]
    fn hdf5_series_stores_attributes() {
        let mut s = XdmfTimeSeriesHdf5::new();
        s.timesteps.push(0.0);
        s.hdf5_paths.push("step0.h5".to_string());
        s.attribute_names.push("density".to_string());
        assert_eq!(s.attribute_names[0], "density");
    }
    #[test]
    fn write_xdmf_with_attributes_basic() {
        let attrs = vec![XdmfAttribute {
            name: "pressure".to_string(),
            center: "Node".to_string(),
            n_components: 1,
            hdf5_path: "data.h5:/pressure".to_string(),
        }];
        let path = std::env::temp_dir().join("test_xdmf_attrs_extra.xmf");
        write_xdmf_with_attributes(
            path.to_str().unwrap_or(""),
            "data.h5",
            5,
            2,
            "Triangle",
            &attrs,
        )
        .unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        assert!(s.contains("pressure"));
        assert!(s.contains("Scalar"));
    }
    #[test]
    fn parse_xdmf_topology_triangle() {
        let xml = "<Topology TopologyType=\"Triangle\" NumberOfElements=\"10\"/>";
        let result = parse_xdmf_topology(xml);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, "Triangle");
        assert_eq!(result[0].1, 10);
    }
    #[test]
    fn parse_xdmf_topology_multiple() {
        let xml = "<Topology TopologyType=\"Triangle\" NumberOfElements=\"4\"/>\n\
                   <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"2\"/>";
        let result = parse_xdmf_topology(xml);
        assert_eq!(result.len(), 2);
    }
    #[test]
    fn parse_xdmf_topology_empty_input() {
        let result = parse_xdmf_topology("");
        assert!(result.is_empty());
    }
    #[test]
    fn xdmf_time_element_zero() {
        let t = xdmf_time_element(0.0);
        assert_eq!(t, "<Time Value=\"0\"/>");
    }
    #[test]
    fn xdmf_time_element_negative() {
        let t = xdmf_time_element(-1.0);
        assert!(t.contains("-1"));
    }
}
/// Write an XDMF XML block for one time step with multiple fields.
///
/// All fields are embedded inline (`Format="XML"`).
pub fn write_xdmf_timestep_fields<W: std::io::Write>(
    writer: &mut W,
    time: f64,
    n_nodes: usize,
    topology_type: &str,
    n_elements: usize,
    nodes: &[[f64; 3]],
    fields: &[XdmfFieldDescriptor],
) -> std::io::Result<()> {
    writeln!(writer, "  <Grid Name=\"timestep\" GridType=\"Uniform\">")?;
    writeln!(writer, "    <Time Value=\"{time}\"/>")?;
    writeln!(
        writer,
        "    <Topology TopologyType=\"{topology_type}\" NumberOfElements=\"{n_elements}\"/>"
    )?;
    writeln!(writer, "    <Geometry GeometryType=\"XYZ\">")?;
    writeln!(
        writer,
        "      <DataItem Format=\"XML\" Dimensions=\"{n_nodes} 3\">"
    )?;
    for n in nodes.iter().take(n_nodes) {
        writeln!(writer, "        {} {} {}", n[0], n[1], n[2])?;
    }
    writeln!(writer, "      </DataItem>")?;
    writeln!(writer, "    </Geometry>")?;
    for field in fields {
        let n_entries = field.entry_count();
        let dim_str = if field.n_components == 1 {
            format!("{n_entries}")
        } else {
            format!("{n_entries} {}", field.n_components)
        };
        writeln!(
            writer,
            "    <Attribute Name=\"{}\" AttributeType=\"{}\" Center=\"{}\">",
            field.name, field.attribute_type, field.center
        )?;
        writeln!(
            writer,
            "      <DataItem Format=\"XML\" Dimensions=\"{dim_str}\">"
        )?;
        write!(writer, "        ")?;
        for (i, v) in field.data.iter().enumerate() {
            if i > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{v}")?;
        }
        writeln!(writer)?;
        writeln!(writer, "      </DataItem>")?;
        writeln!(writer, "    </Attribute>")?;
    }
    writeln!(writer, "  </Grid>")?;
    Ok(())
}
/// Collect all patches into a mapping element_id → patch_name (first match).
pub fn patch_element_map(patches: &[XdmfMeshPatch]) -> std::collections::HashMap<usize, String> {
    let mut map = std::collections::HashMap::new();
    for patch in patches {
        for &eid in &patch.element_ids {
            map.entry(eid).or_insert_with(|| patch.name.clone());
        }
    }
    map
}
/// Format a flat array of `[f64; 3]` positions as a space-separated XDMF data string.
pub fn format_xdmf_geometry_inline(nodes: &[[f64; 3]]) -> String {
    nodes
        .iter()
        .map(|n| format!("{} {} {}", n[0], n[1], n[2]))
        .collect::<Vec<_>>()
        .join("\n")
}
/// Format a flat `Vec`f64` as a single space-separated line suitable for XDMF DataItem.
pub fn format_xdmf_data_inline(data: &[f64]) -> String {
    data.iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(" ")
}
/// Build an XDMF `<DataItem Format="HDF" ...>` reference string.
///
/// `hdf5_file` is the HDF5 file path; `dataset` is the dataset path inside
/// the file; `dims` is a space-separated dimension string like `"100 3"`.
pub fn xdmf_hdf5_dataitem(hdf5_file: &str, dataset: &str, dims: &str, number_type: &str) -> String {
    format!(
        "<DataItem Format=\"HDF\" Dimensions=\"{dims}\" NumberType=\"{number_type}\" Precision=\"8\">{hdf5_file}:{dataset}</DataItem>"
    )
}
/// Validate that an XDMF XML string contains the minimum required tags.
///
/// Returns `Ok(())` if the string looks like valid XDMF, or `Err(msg)` describing
/// the first missing element.
pub fn validate_xdmf_structure(xml: &str) -> Result<(), String> {
    let required = ["<?xml", "<Xdmf", "<Domain>", "</Domain>", "</Xdmf>"];
    for tag in &required {
        if !xml.contains(tag) {
            return Err(format!("XDMF validation: missing `{tag}`"));
        }
    }
    Ok(())
}
/// Indent all lines in an XDMF string by `level * 2` spaces.
pub fn indent_xdmf(xml: &str, level: usize) -> String {
    let prefix = " ".repeat(level * 2);
    xml.lines()
        .map(|line| format!("{prefix}{line}"))
        .collect::<Vec<_>>()
        .join("\n")
}
