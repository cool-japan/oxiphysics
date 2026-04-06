//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_div_ceil, clippy::too_many_arguments)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use std::io::Write;

/// Topology types supported for unstructured mesh output.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum XdmfTopologyType {
    /// Linear triangle (3 nodes).
    Triangle,
    /// Linear tetrahedron (4 nodes).
    Tetrahedron,
    /// Linear hexahedron (8 nodes).
    Hexahedron,
    /// Linear quadrilateral (4 nodes).
    Quadrilateral,
    /// Mixed element types (XDMF Mixed topology).
    Mixed,
}
impl XdmfTopologyType {
    /// Return the XDMF topology type string.
    pub fn xdmf_name(self) -> &'static str {
        match self {
            XdmfTopologyType::Triangle => "Triangle",
            XdmfTopologyType::Tetrahedron => "Tetrahedron",
            XdmfTopologyType::Hexahedron => "Hexahedron",
            XdmfTopologyType::Quadrilateral => "Quadrilateral",
            XdmfTopologyType::Mixed => "Mixed",
        }
    }
    /// Nodes per element (0 for Mixed, since element size varies).
    pub fn nodes_per_element(self) -> usize {
        match self {
            XdmfTopologyType::Triangle => 3,
            XdmfTopologyType::Tetrahedron => 4,
            XdmfTopologyType::Hexahedron => 8,
            XdmfTopologyType::Quadrilateral => 4,
            XdmfTopologyType::Mixed => 0,
        }
    }
}
/// A single time step in an XDMF time series.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfStep {
    /// Simulation time.
    pub time: f64,
    /// Number of points.
    pub n_points: usize,
    /// Position data.
    pub position_data: Vec<[f64; 3]>,
    /// Named scalar fields.
    pub scalar_fields: Vec<(String, Vec<f64>)>,
}
/// Parameters for a structured uniform 3-D grid.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfUniformGrid {
    /// Grid name.
    pub name: String,
    /// Number of nodes along each axis `[nx, ny, nz]`.
    pub dimensions: [usize; 3],
    /// Origin of the grid `[ox, oy, oz]`.
    pub origin: [f64; 3],
    /// Grid spacing `[dx, dy, dz]`.
    pub spacing: [f64; 3],
}
/// Build an XDMF multi-block (spatial collection) document containing
/// multiple uniform grids sharing the same domain.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct XdmfMultiBlock {
    /// Named blocks (each block is a complete XDMF Uniform grid XML fragment,
    /// without the outer `<?xml …>`, `<Xdmf …>`, or ``Domain` wrappers).
    pub(super) blocks: Vec<(String, String)>,
}
#[allow(dead_code)]
impl XdmfMultiBlock {
    /// Create an empty multi-block container.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a named block given its inner grid XML (starting with `<Grid …>`).
    pub fn add_block(&mut self, name: &str, grid_xml: &str) {
        self.blocks.push((name.to_string(), grid_xml.to_string()));
    }
    /// Number of blocks.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }
    /// Returns `true` if no blocks have been added.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }
    /// Serialize to a complete XDMF document.
    pub fn to_xml(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<Xdmf Version=\"3.0\">\n");
        s.push_str("  <Domain>\n");
        s.push_str(
            "    <Grid Name=\"MultiBlock\" GridType=\"Collection\" CollectionType=\"Spatial\">\n",
        );
        for (_name, grid) in &self.blocks {
            for line in grid.lines() {
                s.push_str("      ");
                s.push_str(line);
                s.push('\n');
            }
        }
        s.push_str("    </Grid>\n");
        s.push_str("  </Domain>\n");
        s.push_str("</Xdmf>\n");
        s
    }
}
/// A time series referencing HDF5 datasets by path.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfTimeSeriesHdf5 {
    /// Time values for each step.
    pub timesteps: Vec<f64>,
    /// HDF5 file paths for each step's geometry data.
    pub hdf5_paths: Vec<String>,
    /// Attribute names to reference in HDF5.
    pub attribute_names: Vec<String>,
}
#[allow(dead_code)]
impl XdmfTimeSeriesHdf5 {
    /// Create a new empty HDF5-backed time series.
    pub fn new() -> Self {
        Self {
            timesteps: Vec::new(),
            hdf5_paths: Vec::new(),
            attribute_names: Vec::new(),
        }
    }
    /// Write a multi-timestep XDMF collection to `path`.
    ///
    /// `n_nodes` and `n_elements` describe the mesh size for each step;
    /// `topology` is the XDMF topology type string (e.g., `"Triangle"`).
    pub fn write_collection(
        &self,
        path: &str,
        n_nodes: usize,
        n_elements: usize,
        topology: &str,
    ) -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        writeln!(f, "<?xml version=\"1.0\"?>")?;
        writeln!(f, "<Xdmf Version=\"3.0\">")?;
        writeln!(f, "  <Domain>")?;
        writeln!(
            f,
            "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">"
        )?;
        for (step_idx, &t) in self.timesteps.iter().enumerate() {
            let h5path = self
                .hdf5_paths
                .get(step_idx)
                .map(|s| s.as_str())
                .unwrap_or("data.h5");
            writeln!(
                f,
                "      <Grid Name=\"step{}\" GridType=\"Uniform\">",
                step_idx
            )?;
            writeln!(f, "        <Time Value=\"{}\"/>", t)?;
            writeln!(
                f,
                "        <Topology TopologyType=\"{}\" NumberOfElements=\"{}\"/>",
                topology, n_elements
            )?;
            writeln!(f, "        <Geometry GeometryType=\"XYZ\">")?;
            writeln!(
                f,
                "          <DataItem Format=\"HDF\" Dimensions=\"{} 3\">{}:/coordinates</DataItem>",
                n_nodes, h5path
            )?;
            writeln!(f, "        </Geometry>")?;
            for attr in &self.attribute_names {
                writeln!(
                    f,
                    "        <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">",
                    attr
                )?;
                writeln!(
                    f,
                    "          <DataItem Format=\"HDF\" Dimensions=\"{}\">{}/{}</DataItem>",
                    n_nodes, h5path, attr
                )?;
                writeln!(f, "        </Attribute>")?;
            }
            writeln!(f, "      </Grid>")?;
        }
        writeln!(f, "    </Grid>")?;
        writeln!(f, "  </Domain>")?;
        writeln!(f, "</Xdmf>")?;
        Ok(())
    }
}
/// A structured (regular) Cartesian grid with optional scalar/vector fields.
///
/// This is the primary building block for CFD output on structured meshes.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfStructuredGrid {
    /// Grid name shown in post-processors.
    pub name: String,
    /// Number of nodes along X (nodes = cells + 1 per axis).
    pub ni: usize,
    /// Number of nodes along Y.
    pub nj: usize,
    /// Number of nodes along Z.
    pub nk: usize,
    /// Physical origin of the grid.
    pub origin: [f64; 3],
    /// Grid spacing along X.
    pub dx: f64,
    /// Grid spacing along Y.
    pub dy: f64,
    /// Grid spacing along Z.
    pub dz: f64,
    /// Node-centered scalar fields `(name, flat values, len = ni*nj*nk)`.
    pub node_scalars: Vec<(String, Vec<f64>)>,
    /// Node-centered vector fields `(name, flat values, len = ni*nj*nk)`.
    pub node_vectors: Vec<(String, Vec<[f64; 3]>)>,
    /// Cell-centered scalar fields `(name, flat values, len = (ni-1)*(nj-1)*(nk-1))`.
    pub cell_scalars: Vec<(String, Vec<f64>)>,
}
#[allow(dead_code)]
impl XdmfStructuredGrid {
    /// Create a new structured grid with no fields.
    pub fn new(
        name: &str,
        ni: usize,
        nj: usize,
        nk: usize,
        origin: [f64; 3],
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> Self {
        Self {
            name: name.to_string(),
            ni,
            nj,
            nk,
            origin,
            dx,
            dy,
            dz,
            node_scalars: Vec::new(),
            node_vectors: Vec::new(),
            cell_scalars: Vec::new(),
        }
    }
    /// Total node count.
    pub fn n_nodes(&self) -> usize {
        self.ni * self.nj * self.nk
    }
    /// Total cell count.
    pub fn n_cells(&self) -> usize {
        (self.ni.saturating_sub(1)) * (self.nj.saturating_sub(1)) * (self.nk.saturating_sub(1))
    }
    /// Add a node-centered scalar field.
    pub fn add_node_scalar(&mut self, name: &str, data: Vec<f64>) {
        self.node_scalars.push((name.to_string(), data));
    }
    /// Add a node-centered vector field.
    pub fn add_node_vector(&mut self, name: &str, data: Vec<[f64; 3]>) {
        self.node_vectors.push((name.to_string(), data));
    }
    /// Add a cell-centered scalar field.
    pub fn add_cell_scalar(&mut self, name: &str, data: Vec<f64>) {
        self.cell_scalars.push((name.to_string(), data));
    }
    /// Serialize the structured grid to XDMF XML using ORIGIN_DXDYDZ geometry.
    pub fn to_xml(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<Xdmf Version=\"3.0\">\n");
        s.push_str("  <Domain>\n");
        s.push_str(&format!(
            "    <Grid Name=\"{}\" GridType=\"Uniform\">\n",
            self.name
        ));
        s.push_str(&format!(
            "      <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"{} {} {}\"/>\n",
            self.nk, self.nj, self.ni
        ));
        s.push_str("      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
        s.push_str(&format!(
            "        <DataItem Format=\"XML\" Dimensions=\"3\">{} {} {}</DataItem>\n",
            self.origin[2], self.origin[1], self.origin[0]
        ));
        s.push_str(&format!(
            "        <DataItem Format=\"XML\" Dimensions=\"3\">{} {} {}</DataItem>\n",
            self.dz, self.dy, self.dx
        ));
        s.push_str("      </Geometry>\n");
        for (name, data) in &self.node_scalars {
            s.push_str(&format!(
                "      <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">\n",
                name
            ));
            s.push_str(&format!(
                "        <DataItem Format=\"XML\" Dimensions=\"{}\">\n          ",
                data.len()
            ));
            for (i, v) in data.iter().enumerate() {
                if i > 0 {
                    s.push(' ');
                }
                s.push_str(&format!("{}", v));
            }
            s.push_str("\n        </DataItem>\n      </Attribute>\n");
        }
        for (name, data) in &self.node_vectors {
            s.push_str(&format!(
                "      <Attribute Name=\"{}\" AttributeType=\"Vector\" Center=\"Node\">\n",
                name
            ));
            s.push_str(&format!(
                "        <DataItem Format=\"XML\" Dimensions=\"{} 3\">\n",
                data.len()
            ));
            for v in data {
                s.push_str(&format!("          {} {} {}\n", v[0], v[1], v[2]));
            }
            s.push_str("        </DataItem>\n      </Attribute>\n");
        }
        for (name, data) in &self.cell_scalars {
            s.push_str(&format!(
                "      <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
                name
            ));
            s.push_str(&format!(
                "        <DataItem Format=\"XML\" Dimensions=\"{}\">\n          ",
                data.len()
            ));
            for (i, v) in data.iter().enumerate() {
                if i > 0 {
                    s.push(' ');
                }
                s.push_str(&format!("{}", v));
            }
            s.push_str("\n        </DataItem>\n      </Attribute>\n");
        }
        s.push_str("    </Grid>\n");
        s.push_str("  </Domain>\n");
        s.push_str("</Xdmf>\n");
        s
    }
}
/// A stateful writer that accumulates XDMF XML into a `String`.
///
/// Supports opening/closing domains, grids, topology, geometry, and attributes
/// incrementally, allowing complex documents to be built programmatically.
#[allow(dead_code)]
pub struct XdmfWriter {
    pub(super) buf: String,
    pub(super) indent: usize,
}
#[allow(dead_code)]
impl XdmfWriter {
    /// Create a new writer and emit the XML declaration and root ``Xdmf` element.
    pub fn new() -> Self {
        let mut w = XdmfWriter {
            buf: String::new(),
            indent: 0,
        };
        w.buf.push_str("<?xml version=\"1.0\"?>\n");
        w.buf.push_str("<Xdmf Version=\"3.0\">\n");
        w.indent = 1;
        w
    }
    fn push_line(&mut self, s: &str) {
        let pad = "  ".repeat(self.indent);
        self.buf.push_str(&pad);
        self.buf.push_str(s);
        self.buf.push('\n');
    }
    /// Open ``Domain`.
    pub fn open_domain(&mut self) {
        self.push_line("<Domain>");
        self.indent += 1;
    }
    /// Close `</Domain>`.
    pub fn close_domain(&mut self) {
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Domain>");
    }
    /// Open a ``Grid` element with `name` and optional `GridType`.
    pub fn open_grid(&mut self, name: &str, grid_type: &str) {
        self.push_line(&format!(
            "<Grid Name=\"{}\" GridType=\"{}\">",
            name, grid_type
        ));
        self.indent += 1;
    }
    /// Open a temporal collection grid.
    pub fn open_temporal_collection(&mut self, name: &str) {
        self.push_line(&format!(
            "<Grid Name=\"{}\" GridType=\"Collection\" CollectionType=\"Temporal\">",
            name
        ));
        self.indent += 1;
    }
    /// Emit a `<Time Value="t"/>` element.
    pub fn write_time(&mut self, t: f64) {
        self.push_line(&format!("<Time Value=\"{}\"/>", t));
    }
    /// Close `</Grid>`.
    pub fn close_grid(&mut self) {
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Grid>");
    }
    /// Emit a ``Topology` element for a polyvertex (point cloud).
    pub fn write_polyvertex_topology(&mut self, n: usize) {
        self.push_line(&format!(
            "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"{}\"/>",
            n
        ));
    }
    /// Emit a ``Topology` element for an unstructured mesh.
    pub fn write_unstructured_topology(
        &mut self,
        topo_type: &str,
        n_elements: usize,
        connectivity: &[usize],
        npe: usize,
    ) {
        self.push_line(&format!(
            "<Topology TopologyType=\"{}\" NumberOfElements=\"{}\">",
            topo_type, n_elements
        ));
        self.indent += 1;
        self.push_line(&format!(
            "<DataItem Format=\"XML\" Dimensions=\"{} {}\">\n{}  ",
            n_elements,
            npe,
            "  ".repeat(self.indent)
        ));
        let row_strings: Vec<String> = connectivity
            .chunks(npe)
            .map(|chunk| {
                chunk
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect();
        let pad = "  ".repeat(self.indent);
        for r in &row_strings {
            self.buf.push_str(&pad);
            self.buf.push_str(r);
            self.buf.push('\n');
        }
        self.push_line("</DataItem>");
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Topology>");
    }
    /// Emit a `<Geometry GeometryType="XYZ">` with inline node coordinates.
    pub fn write_xyz_geometry(&mut self, nodes: &[[f64; 3]]) {
        self.push_line("<Geometry GeometryType=\"XYZ\">");
        self.indent += 1;
        self.push_line(&format!(
            "<DataItem Format=\"XML\" Dimensions=\"{} 3\">",
            nodes.len()
        ));
        self.indent += 1;
        let pad = "  ".repeat(self.indent);
        for p in nodes {
            self.buf.push_str(&pad);
            self.buf.push_str(&format!("{} {} {}\n", p[0], p[1], p[2]));
        }
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</DataItem>");
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Geometry>");
    }
    /// Emit a scalar attribute with inline data.
    pub fn write_scalar_attribute(&mut self, name: &str, center: &str, values: &[f64]) {
        self.push_line(&format!(
            "<Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"{}\">",
            name, center
        ));
        self.indent += 1;
        self.push_line(&format!(
            "<DataItem Format=\"XML\" Dimensions=\"{}\">",
            values.len()
        ));
        self.indent += 1;
        let pad = "  ".repeat(self.indent);
        self.buf.push_str(&pad);
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                self.buf.push(' ');
            }
            self.buf.push_str(&format!("{}", v));
        }
        self.buf.push('\n');
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</DataItem>");
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Attribute>");
    }
    /// Emit a vector attribute with inline data.
    pub fn write_vector_attribute(&mut self, name: &str, center: &str, vectors: &[[f64; 3]]) {
        self.push_line(&format!(
            "<Attribute Name=\"{}\" AttributeType=\"Vector\" Center=\"{}\">",
            name, center
        ));
        self.indent += 1;
        self.push_line(&format!(
            "<DataItem Format=\"XML\" Dimensions=\"{} 3\">",
            vectors.len()
        ));
        self.indent += 1;
        let pad = "  ".repeat(self.indent);
        for v in vectors {
            self.buf.push_str(&pad);
            self.buf.push_str(&format!("{} {} {}\n", v[0], v[1], v[2]));
        }
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</DataItem>");
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Attribute>");
    }
    /// Emit a HDF5-reference attribute.
    pub fn write_hdf5_attribute(
        &mut self,
        name: &str,
        center: &str,
        attr_type: &str,
        dims: &str,
        hdf5_ref: &str,
    ) {
        self.push_line(&format!(
            "<Attribute Name=\"{}\" AttributeType=\"{}\" Center=\"{}\">",
            name, attr_type, center
        ));
        self.indent += 1;
        self.push_line(&format!(
            "<DataItem Format=\"HDF\" Dimensions=\"{}\">{}</DataItem>",
            dims, hdf5_ref
        ));
        self.indent = self.indent.saturating_sub(1);
        self.push_line("</Attribute>");
    }
    /// Emit the closing `</Xdmf>` element and return the accumulated XML string.
    pub fn finish(mut self) -> String {
        self.buf.push_str("</Xdmf>\n");
        self.buf
    }
    /// Return the current buffer without closing (useful for inspection).
    pub fn peek(&self) -> &str {
        &self.buf
    }
}
/// Builds a canonical HDF5 dataset path for use in XDMF ``DataItem` elements.
///
/// A path looks like `file.h5:/group/dataset`.
#[allow(dead_code)]
pub struct Hdf5DataItemBuilder {
    pub(super) filename: String,
    pub(super) group: String,
}
#[allow(dead_code)]
impl Hdf5DataItemBuilder {
    /// Create a builder targeting `filename`.
    pub fn new(filename: &str) -> Self {
        Self {
            filename: filename.to_owned(),
            group: "/".to_owned(),
        }
    }
    /// Set the HDF5 group prefix (must start with `/`).
    pub fn group(mut self, group: &str) -> Self {
        self.group = group.to_owned();
        self
    }
    /// Build a ``DataItem` XML fragment referencing `dataset_name` under
    /// the configured group, with the given `dimensions` string (e.g., `"100 3"`).
    pub fn build(&self, dataset_name: &str, dimensions: &str) -> String {
        let path = if self.group == "/" {
            format!("/{}", dataset_name)
        } else {
            format!("{}/{}", self.group, dataset_name)
        };
        format!(
            "<DataItem Format=\"HDF\" Dimensions=\"{}\">\n  {}:{}\n</DataItem>",
            dimensions, self.filename, path
        )
    }
}
/// A named grid inside a domain collection.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfGridEntry {
    /// Grid name.
    pub name: String,
    /// Time value for this grid.
    pub time: f64,
    /// Node positions.
    pub nodes: Vec<[f64; 3]>,
    /// Flat connectivity array.
    pub connectivity: Vec<usize>,
    /// Topology type string.
    pub topology_type: String,
}
impl XdmfGridEntry {
    /// Number of nodes in this grid.
    #[allow(dead_code)]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Compute the axis-aligned bounding box of all nodes.
    /// Returns `([min_x, min_y, min_z], [max_x, max_y, max_z])`.
    #[allow(dead_code)]
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for n in &self.nodes {
            for k in 0..3 {
                if n[k] < lo[k] {
                    lo[k] = n[k];
                }
                if n[k] > hi[k] {
                    hi[k] = n[k];
                }
            }
        }
        (lo, hi)
    }
    /// Compute the centroid of all nodes.
    #[allow(dead_code)]
    pub fn centroid(&self) -> [f64; 3] {
        if self.nodes.is_empty() {
            return [0.0; 3];
        }
        let mut s = [0.0_f64; 3];
        for n in &self.nodes {
            s[0] += n[0];
            s[1] += n[1];
            s[2] += n[2];
        }
        let inv = 1.0 / self.nodes.len() as f64;
        [s[0] * inv, s[1] * inv, s[2] * inv]
    }
}
/// Describes an XDMF attribute (field) referencing HDF5 data.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfAttribute {
    /// Attribute name.
    pub name: String,
    /// Centering: `"Node"`, `"Cell"`, etc.
    pub center: String,
    /// Number of components (1 = scalar, 3 = vector, etc.).
    pub n_components: usize,
    /// HDF5 dataset path (e.g., `"data.h5:/temperature"`).
    pub hdf5_path: String,
}
/// A single time step for an unstructured mesh time series.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct XdmfMeshStep {
    /// Simulation time.
    pub time: f64,
    /// Node positions (3-D).
    pub nodes: Vec<[f64; 3]>,
    /// Flat connectivity array (length = n_elements * nodes_per_element).
    pub connectivity: Vec<usize>,
    /// Element topology.
    pub topology: XdmfTopologyType,
    /// Named node-centered scalar fields.
    pub node_scalars: Vec<(String, Vec<f64>)>,
    /// Named node-centered vector fields (each vec is a `[f64;3]` per node).
    pub node_vectors: Vec<(String, Vec<[f64; 3]>)>,
}
impl XdmfMeshStep {
    /// Number of elements in this step.
    pub fn n_elements(&self) -> usize {
        let npe = self.topology.nodes_per_element();
        self.connectivity.len().checked_div(npe).unwrap_or(0)
    }
}
/// A collection of grids forming an XDMF domain.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct XdmfDomainCollection {
    /// All grids in this collection.
    pub grids: Vec<XdmfGridEntry>,
}
impl XdmfDomainCollection {
    /// Create an empty domain collection.
    #[allow(dead_code)]
    pub fn new() -> Self {
        XdmfDomainCollection { grids: Vec::new() }
    }
    /// Add a grid entry to this collection.
    #[allow(dead_code)]
    pub fn add_grid(&mut self, grid: XdmfGridEntry) {
        self.grids.push(grid);
    }
    /// Total number of nodes across all grids.
    #[allow(dead_code)]
    pub fn total_node_count(&self) -> usize {
        self.grids.iter().map(|g| g.node_count()).sum()
    }
    /// Find a grid by name.
    #[allow(dead_code)]
    pub fn find_grid(&self, name: &str) -> Option<&XdmfGridEntry> {
        self.grids.iter().find(|g| g.name == name)
    }
    /// Write the full XDMF XML for this domain collection.
    #[allow(dead_code)]
    pub fn write_xml<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(writer, "<?xml version=\"1.0\"?>")?;
        writeln!(writer, "<Xdmf Version=\"3.0\">")?;
        writeln!(writer, "  <Domain>")?;
        writeln!(
            writer,
            "    <Grid Name=\"collection\" GridType=\"Collection\" CollectionType=\"Temporal\">"
        )?;
        for grid in &self.grids {
            let nn = grid.node_count();
            writeln!(
                writer,
                "      <Grid Name=\"{}\" GridType=\"Uniform\">",
                grid.name
            )?;
            writeln!(writer, "        <Time Value=\"{}\"/>", grid.time)?;
            writeln!(
                writer,
                "        <Topology TopologyType=\"{}\" NumberOfElements=\"{}\"/>",
                grid.topology_type,
                grid.connectivity.len()
            )?;
            writeln!(writer, "        <Geometry GeometryType=\"XYZ\">")?;
            writeln!(
                writer,
                "          <DataItem Format=\"XML\" Dimensions=\"{nn} 3\">"
            )?;
            for n in &grid.nodes {
                writeln!(writer, "            {} {} {}", n[0], n[1], n[2])?;
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
    /// Time range: (min_time, max_time).  Returns (0.0, 0.0) if no grids.
    #[allow(dead_code)]
    pub fn time_range(&self) -> (f64, f64) {
        if self.grids.is_empty() {
            return (0.0, 0.0);
        }
        let t_min = self
            .grids
            .iter()
            .map(|g| g.time)
            .fold(f64::INFINITY, f64::min);
        let t_max = self
            .grids
            .iter()
            .map(|g| g.time)
            .fold(f64::NEG_INFINITY, f64::max);
        (t_min, t_max)
    }
}
/// Describes a single data field to embed inline in XDMF XML.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfFieldDescriptor {
    /// Name of the field (e.g. `"pressure"`, `"velocity"`).
    pub name: String,
    /// Attribute type: `"Scalar"`, `"Vector"`, or `"Tensor"`.
    pub attribute_type: String,
    /// Center: `"Node"`, `"Cell"`, or `"Grid"`.
    pub center: String,
    /// Flat data values (row-major).
    pub data: Vec<f64>,
    /// Number of components per entry (1 for scalar, 3 for vector, 9 for tensor).
    pub n_components: usize,
}
impl XdmfFieldDescriptor {
    /// Create a new scalar field descriptor.
    #[allow(dead_code)]
    pub fn scalar(name: &str, data: Vec<f64>) -> Self {
        XdmfFieldDescriptor {
            name: name.to_string(),
            attribute_type: "Scalar".to_string(),
            center: "Node".to_string(),
            data,
            n_components: 1,
        }
    }
    /// Create a new vector field descriptor (3 components per node).
    #[allow(dead_code)]
    pub fn vector(name: &str, data: Vec<f64>) -> Self {
        XdmfFieldDescriptor {
            name: name.to_string(),
            attribute_type: "Vector".to_string(),
            center: "Node".to_string(),
            data,
            n_components: 3,
        }
    }
    /// Number of logical entries (nodes/cells) this field covers.
    #[allow(dead_code)]
    pub fn entry_count(&self) -> usize {
        self.data.len().checked_div(self.n_components).unwrap_or(0)
    }
    /// Compute the Lp-norm of all data values.
    #[allow(dead_code)]
    pub fn data_lp_norm(&self, p: f64) -> f64 {
        if p <= 0.0 || self.data.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.data.iter().map(|v| v.abs().powf(p)).sum();
        sum.powf(1.0 / p)
    }
    /// Return the maximum absolute value in the field data.
    #[allow(dead_code)]
    pub fn max_abs(&self) -> f64 {
        self.data
            .iter()
            .cloned()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()))
    }
    /// Return the minimum value in the field data (or 0 if empty).
    #[allow(dead_code)]
    pub fn min_value(&self) -> f64 {
        self.data.iter().cloned().fold(f64::INFINITY, f64::min)
    }
    /// Return the maximum value in the field data (or 0 if empty).
    #[allow(dead_code)]
    pub fn max_value(&self) -> f64 {
        self.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
}
/// A time-varying unstructured mesh collection (multi-topology supported across
/// steps, though a single topology per step is assumed).
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct XdmfMeshTimeSeries {
    /// Ordered list of mesh steps.
    pub steps: Vec<XdmfMeshStep>,
}
#[allow(dead_code)]
impl XdmfMeshTimeSeries {
    /// Create an empty mesh time series.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append one mesh time step.
    pub fn add_step(&mut self, step: XdmfMeshStep) {
        self.steps.push(step);
    }
    /// Number of steps.
    pub fn len(&self) -> usize {
        self.steps.len()
    }
    /// Returns `true` if no steps have been added.
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }
    /// Return the simulation time values of all steps.
    pub fn times(&self) -> Vec<f64> {
        self.steps.iter().map(|s| s.time).collect()
    }
    /// Serialize the mesh time series to XDMF XML.
    pub fn to_xml(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<Xdmf Version=\"3.0\">\n");
        s.push_str("  <Domain>\n");
        s.push_str(
            "    <Grid Name=\"MeshTimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n",
        );
        for step in &self.steps {
            let n_nodes = step.nodes.len();
            let n_elements = step.n_elements();
            let npe = step.topology.nodes_per_element();
            s.push_str("      <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
            s.push_str(&format!("        <Time Value=\"{}\"/>\n", step.time));
            s.push_str(&format!(
                "        <Topology TopologyType=\"{}\" NumberOfElements=\"{}\">\n",
                step.topology.xdmf_name(),
                n_elements
            ));
            s.push_str(&format!(
                "          <DataItem Format=\"XML\" Dimensions=\"{} {}\">\n",
                n_elements, npe
            ));
            for chunk in step.connectivity.chunks(npe) {
                let row: Vec<String> = chunk.iter().map(|&i| i.to_string()).collect();
                s.push_str(&format!("            {}\n", row.join(" ")));
            }
            s.push_str("          </DataItem>\n");
            s.push_str("        </Topology>\n");
            s.push_str("        <Geometry GeometryType=\"XYZ\">\n");
            s.push_str(&format!(
                "          <DataItem Format=\"XML\" Dimensions=\"{} 3\">\n",
                n_nodes
            ));
            for p in &step.nodes {
                s.push_str(&format!("            {} {} {}\n", p[0], p[1], p[2]));
            }
            s.push_str("          </DataItem>\n");
            s.push_str("        </Geometry>\n");
            for (name, values) in &step.node_scalars {
                s.push_str(&format!(
                    "        <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">\n",
                    name
                ));
                s.push_str(&format!(
                    "          <DataItem Format=\"XML\" Dimensions=\"{}\">\n",
                    values.len()
                ));
                s.push_str("            ");
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        s.push(' ');
                    }
                    s.push_str(&format!("{}", v));
                }
                s.push('\n');
                s.push_str("          </DataItem>\n");
                s.push_str("        </Attribute>\n");
            }
            for (name, vdata) in &step.node_vectors {
                s.push_str(&format!(
                    "        <Attribute Name=\"{}\" AttributeType=\"Vector\" Center=\"Node\">\n",
                    name
                ));
                s.push_str(&format!(
                    "          <DataItem Format=\"XML\" Dimensions=\"{} 3\">\n",
                    vdata.len()
                ));
                for v in vdata {
                    s.push_str(&format!("            {} {} {}\n", v[0], v[1], v[2]));
                }
                s.push_str("          </DataItem>\n");
                s.push_str("        </Attribute>\n");
            }
            s.push_str("      </Grid>\n");
        }
        s.push_str("    </Grid>\n");
        s.push_str("  </Domain>\n");
        s.push_str("</Xdmf>\n");
        s
    }
}
/// An XDMF time series collection.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfTimeSeries {
    /// Time steps.
    pub steps: Vec<XdmfStep>,
}
#[allow(dead_code)]
impl XdmfTimeSeries {
    /// Create a new empty time series.
    pub fn new() -> Self {
        Self { steps: Vec::new() }
    }
    /// Add a time step with positions and optional scalar fields.
    pub fn add_step(
        &mut self,
        time: f64,
        positions: Vec<[f64; 3]>,
        scalars: Vec<(String, Vec<f64>)>,
    ) {
        let n_points = positions.len();
        self.steps.push(XdmfStep {
            time,
            n_points,
            position_data: positions,
            scalar_fields: scalars,
        });
    }
    /// Serialize the time series to XDMF XML.
    pub fn to_xml(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<Xdmf Version=\"3.0\">\n");
        s.push_str("  <Domain>\n");
        s.push_str(
            "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n",
        );
        for step in &self.steps {
            s.push_str("      <Grid Name=\"particles\" GridType=\"Uniform\">\n");
            s.push_str(&format!("        <Time Value=\"{}\"/>\n", step.time));
            s.push_str(&format!(
                "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"{}\"/>\n",
                step.n_points
            ));
            s.push_str("        <Geometry GeometryType=\"XYZ\">\n");
            s.push_str(&format!(
                "          <DataItem Format=\"XML\" Dimensions=\"{} 3\">\n",
                step.n_points
            ));
            for p in &step.position_data {
                s.push_str(&format!("            {} {} {}\n", p[0], p[1], p[2]));
            }
            s.push_str("          </DataItem>\n");
            s.push_str("        </Geometry>\n");
            for (name, values) in &step.scalar_fields {
                s.push_str(&format!(
                    "        <Attribute Name=\"{}\" AttributeType=\"Scalar\" Center=\"Node\">\n",
                    name
                ));
                s.push_str(&format!(
                    "          <DataItem Format=\"XML\" Dimensions=\"{}\">\n",
                    values.len()
                ));
                s.push_str("            ");
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        s.push(' ');
                    }
                    s.push_str(&format!("{}", v));
                }
                s.push('\n');
                s.push_str("          </DataItem>\n");
                s.push_str("        </Attribute>\n");
            }
            s.push_str("      </Grid>\n");
        }
        s.push_str("    </Grid>\n");
        s.push_str("  </Domain>\n");
        s.push_str("</Xdmf>\n");
        s
    }
}
impl XdmfTimeSeries {
    /// Add a step that includes both scalar and vector fields.
    ///
    /// `vectors` is a list of `(name, data)` pairs where each data entry is
    /// a `[f64; 3]` per node.
    #[allow(dead_code)]
    pub fn add_step_with_vectors(
        &mut self,
        time: f64,
        positions: Vec<[f64; 3]>,
        scalars: Vec<(String, Vec<f64>)>,
        vectors: Vec<(String, Vec<[f64; 3]>)>,
    ) {
        let n_points = positions.len();
        let mut all_scalars = scalars;
        for (vname, vdata) in vectors {
            let mut xs = Vec::with_capacity(n_points);
            let mut ys = Vec::with_capacity(n_points);
            let mut zs = Vec::with_capacity(n_points);
            for v in &vdata {
                xs.push(v[0]);
                ys.push(v[1]);
                zs.push(v[2]);
            }
            all_scalars.push((format!("{}_x", vname), xs));
            all_scalars.push((format!("{}_y", vname), ys));
            all_scalars.push((format!("{}_z", vname), zs));
        }
        self.steps.push(XdmfStep {
            time,
            n_points,
            position_data: positions,
            scalar_fields: all_scalars,
        });
    }
    /// Return the time values of all steps.
    #[allow(dead_code)]
    pub fn times(&self) -> Vec<f64> {
        self.steps.iter().map(|s| s.time).collect()
    }
    /// Return the total number of particles across all steps (sum).
    #[allow(dead_code)]
    pub fn total_particle_count(&self) -> usize {
        self.steps.iter().map(|s| s.n_points).sum()
    }
}
impl XdmfTimeSeries {
    /// Add a frame with positions and no scalar fields.
    ///
    /// This is a convenience alias for [`add_step`](XdmfTimeSeries::add_step)
    /// that uses the common `add_frame` naming convention.
    #[allow(dead_code)]
    pub fn add_frame(&mut self, time: f64, positions: Vec<[f64; 3]>) {
        self.add_step(time, positions, vec![]);
    }
    /// Add a frame with positions and named scalar attributes.
    #[allow(dead_code)]
    pub fn add_frame_with_scalars(
        &mut self,
        time: f64,
        positions: Vec<[f64; 3]>,
        scalars: Vec<(String, Vec<f64>)>,
    ) {
        self.add_step(time, positions, scalars);
    }
    /// Write the time series as XDMF XML to any [`Write`] sink.
    ///
    /// This is the `write_xml` variant of [`to_xml`](XdmfTimeSeries::to_xml)
    /// that streams output rather than building a full string first.
    #[allow(dead_code)]
    pub fn write_xml<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let xml = self.to_xml();
        writer.write_all(xml.as_bytes())
    }
    /// Write the time series to a file at `path`.
    #[allow(dead_code)]
    pub fn write_xml_to_file(&self, path: &str) -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        self.write_xml(&mut f)
    }
}
/// A named patch/region inside a larger mesh, identified by element indices.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct XdmfMeshPatch {
    /// Patch name (e.g. `"inlet"`, `"wall"`).
    pub name: String,
    /// Zero-based element indices belonging to this patch.
    pub element_ids: Vec<usize>,
    /// Optional scalar tag (e.g. boundary condition id).
    pub tag: Option<i32>,
}
impl XdmfMeshPatch {
    /// Create a new patch with a name and element list.
    #[allow(dead_code)]
    pub fn new(name: &str, element_ids: Vec<usize>) -> Self {
        XdmfMeshPatch {
            name: name.to_string(),
            element_ids,
            tag: None,
        }
    }
    /// Number of elements in this patch.
    #[allow(dead_code)]
    pub fn element_count(&self) -> usize {
        self.element_ids.len()
    }
    /// Check if a given element index belongs to this patch.
    #[allow(dead_code)]
    pub fn contains_element(&self, idx: usize) -> bool {
        self.element_ids.contains(&idx)
    }
    /// Merge another patch's elements into this one (duplicates removed).
    #[allow(dead_code)]
    pub fn merge(&mut self, other: &XdmfMeshPatch) {
        for &id in &other.element_ids {
            if !self.element_ids.contains(&id) {
                self.element_ids.push(id);
            }
        }
    }
    /// Generate a CDL-like string representation for debugging.
    #[allow(dead_code)]
    pub fn to_debug_string(&self) -> String {
        format!(
            "patch \"{}\" [{} elements] tag={:?}",
            self.name,
            self.element_ids.len(),
            self.tag
        )
    }
}
/// Basic XDMF XML reader.
#[allow(dead_code)]
pub struct XdmfReader;
#[allow(dead_code)]
impl XdmfReader {
    /// Parse XDMF XML data into a list of time steps.
    ///
    /// This is a simplified parser that looks for `<Time Value="..."/>` and
    /// `NumberOfElements="..."` tags. It does not handle all XDMF features.
    pub fn from_xml(data: &str) -> Result<Vec<XdmfStep>, std::io::Error> {
        let mut steps = Vec::new();
        let mut current_time: Option<f64> = None;
        let mut current_n: Option<usize> = None;
        for line in data.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with("<Time")
                && let Some(start) = trimmed.find("Value=\"")
            {
                let rest = &trimmed[start + 7..];
                if let Some(end) = rest.find('"')
                    && let Ok(t) = rest[..end].parse::<f64>()
                {
                    current_time = Some(t);
                }
            }
            if trimmed.contains("NumberOfElements=\"")
                && let Some(start) = trimmed.find("NumberOfElements=\"")
            {
                let rest = &trimmed[start + 18..];
                if let Some(end) = rest.find('"')
                    && let Ok(n) = rest[..end].parse::<usize>()
                {
                    current_n = Some(n);
                }
            }
            if trimmed == "</Grid>"
                && let (Some(t), Some(n)) = (current_time.take(), current_n.take())
            {
                steps.push(XdmfStep {
                    time: t,
                    n_points: n,
                    position_data: Vec::new(),
                    scalar_fields: Vec::new(),
                });
            }
        }
        Ok(steps)
    }
}
/// A simple XDMF schema descriptor: expected topology type and attribute names.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct XdmfSchema {
    /// Expected topology type string (e.g., `"Polyvertex"`, `"Triangle"`).
    pub topology_type: String,
    /// Attribute names that must be present.
    pub required_attributes: Vec<String>,
}
impl XdmfSchema {
    /// Create a schema.
    pub fn new(topology_type: &str, required_attributes: Vec<String>) -> Self {
        Self {
            topology_type: topology_type.to_owned(),
            required_attributes,
        }
    }
    /// Validate an XDMF XML string against this schema.
    ///
    /// Returns a list of human-readable errors; empty if valid.
    pub fn validate(&self, xml: &str) -> Vec<String> {
        let mut errors = Vec::new();
        if !xml.contains(&format!("TopologyType=\"{}\"", self.topology_type)) {
            errors.push(format!("missing TopologyType=\"{}\"", self.topology_type));
        }
        for attr in &self.required_attributes {
            if !xml.contains(&format!("Name=\"{}\"", attr)) {
                errors.push(format!("missing Attribute Name=\"{}\"", attr));
            }
        }
        errors
    }
}
