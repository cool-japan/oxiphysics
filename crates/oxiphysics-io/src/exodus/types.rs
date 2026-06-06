//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// The primary mesh data structure for Exodus II format.
pub struct ExodusMesh {
    /// Node coordinates as (x, y, z) triples.
    pub coordinates: Vec<[f64; 3]>,
    /// Element blocks.
    pub blocks: Vec<ExodusBlock>,
    /// Node sets.
    pub node_sets: Vec<ExodusNodeSet>,
    /// Side sets.
    pub side_sets: Vec<ExodusSideSet>,
    /// Mesh title string.
    pub title: String,
}
impl ExodusMesh {
    /// Create a new empty mesh with the given title.
    pub fn new(title: &str) -> Self {
        ExodusMesh {
            coordinates: Vec::new(),
            blocks: Vec::new(),
            node_sets: Vec::new(),
            side_sets: Vec::new(),
            title: title.to_string(),
        }
    }
    /// Add a node and return its 1-based ID.
    pub fn add_node(&mut self, x: f64, y: f64, z: f64) -> usize {
        self.coordinates.push([x, y, z]);
        self.coordinates.len()
    }
    /// Add an element block and return its index in `self.blocks`.
    pub fn add_block(
        &mut self,
        id: u32,
        name: &str,
        etype: &str,
        n_nodes_per_elem: usize,
    ) -> usize {
        self.blocks.push(ExodusBlock {
            id,
            name: name.to_string(),
            element_type: etype.to_string(),
            n_nodes_per_elem,
            connectivity: Vec::new(),
        });
        self.blocks.len() - 1
    }
    /// Add an element to the specified block using the provided node IDs (1-based).
    pub fn add_element_to_block(&mut self, block_idx: usize, node_ids: &[usize]) {
        self.blocks[block_idx]
            .connectivity
            .extend_from_slice(node_ids);
    }
    /// Add a node set.
    pub fn add_node_set(&mut self, id: u32, name: &str, nodes: Vec<usize>) {
        self.node_sets.push(ExodusNodeSet {
            id,
            name: name.to_string(),
            node_ids: nodes,
        });
    }
    /// Return the total number of nodes.
    pub fn node_count(&self) -> usize {
        self.coordinates.len()
    }
    /// Return the total number of elements across all blocks.
    pub fn element_count(&self) -> usize {
        self.blocks
            .iter()
            .map(|b| {
                b.connectivity
                    .len()
                    .checked_div(b.n_nodes_per_elem)
                    .unwrap_or(0)
            })
            .sum()
    }
}
/// A complete Exodus result combining mesh topology with time-dependent data.
pub struct ExodusResult {
    /// The mesh topology.
    pub mesh: ExodusMesh,
    /// All time steps with associated variable data.
    pub time_steps: Vec<ExodusTimeStep>,
}
impl ExodusResult {
    /// Generate a human-readable summary string.
    pub fn write_summary(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("Title: {}\n", self.mesh.title));
        s.push_str(&format!("Nodes: {}\n", self.mesh.node_count()));
        s.push_str(&format!("Elements: {}\n", self.mesh.element_count()));
        s.push_str(&format!("Blocks: {}\n", self.mesh.blocks.len()));
        s.push_str(&format!("Node Sets: {}\n", self.mesh.node_sets.len()));
        s.push_str(&format!("Side Sets: {}\n", self.mesh.side_sets.len()));
        s.push_str(&format!("Time Steps: {}\n", self.time_steps.len()));
        for (i, ts) in self.time_steps.iter().enumerate() {
            s.push_str(&format!(
                "  Step {}: t={:.6}, vars={}\n",
                i,
                ts.time,
                ts.variables.len()
            ));
        }
        s
    }
}
/// An element block containing elements of the same type.
pub struct ExodusBlock {
    /// Unique identifier for this block.
    pub id: u32,
    /// Human-readable name.
    pub name: String,
    /// Element type string, e.g. "HEX8", "TET4", "QUAD4".
    pub element_type: String,
    /// Number of nodes per element.
    pub n_nodes_per_elem: usize,
    /// Flattened connectivity array (node IDs, 1-based).
    pub connectivity: Vec<usize>,
}
/// Quality metric result for a single element.
#[derive(Debug, Clone)]
pub struct ElementQuality {
    /// Block id the element belongs to.
    pub block_id: u32,
    /// Zero-based element index within the block.
    pub element_index: usize,
    /// Aspect ratio (min_edge / max_edge for simplices).
    pub aspect_ratio: f64,
    /// Scaled Jacobian for volume elements (1.0 = ideal).
    pub scaled_jacobian: f64,
    /// Minimum edge length.
    pub min_edge: f64,
    /// Maximum edge length.
    pub max_edge: f64,
}
/// Partition result: a list of sub-meshes and their owning partition ids.
pub struct MeshPartition {
    /// Partition index (0-based).
    pub partition_id: usize,
    /// The sub-mesh for this partition.
    pub mesh: ExodusMesh,
}
/// Result of mesh validation.
#[derive(Debug, Clone, Default)]
pub struct MeshValidationReport {
    /// True if no errors were found.
    pub is_valid: bool,
    /// List of error messages.
    pub errors: Vec<String>,
    /// List of warning messages.
    pub warnings: Vec<String>,
}
impl MeshValidationReport {
    pub(super) fn add_error(&mut self, msg: String) {
        self.errors.push(msg);
        self.is_valid = false;
    }
    pub(super) fn add_warning(&mut self, msg: String) {
        self.warnings.push(msg);
    }
}
/// A named set of nodes (boundary condition set).
pub struct ExodusNodeSet {
    /// Unique identifier for this node set.
    pub id: u32,
    /// Human-readable name.
    pub name: String,
    /// List of node IDs (1-based) in this set.
    pub node_ids: Vec<usize>,
}
/// Statistics over a scalar field defined on nodes or elements.
#[derive(Debug, Clone)]
pub struct FieldStats {
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// Population standard deviation.
    pub std_dev: f64,
    /// Number of values.
    pub count: usize,
}
/// A time step containing variable data.
pub struct ExodusTimeStep {
    /// Simulation time for this step.
    pub time: f64,
    /// Variables defined at this time step.
    pub variables: Vec<ExodusVariable>,
}
/// A named set of element sides (face boundary conditions).
pub struct ExodusSideSet {
    /// Unique identifier for this side set.
    pub id: u32,
    /// Human-readable name.
    pub name: String,
    /// Element IDs (1-based) for each side entry.
    pub element_ids: Vec<usize>,
    /// Local side IDs within each element.
    pub side_ids: Vec<usize>,
}
/// A global simulation variable (scalar value valid for the whole domain).
pub struct ExodusGlobalVariable {
    /// Variable name.
    pub name: String,
    /// Simulation time.
    pub time: f64,
    /// Scalar value.
    pub value: f64,
}
/// A single scalar/vector variable attached to an entity (node or element).
pub struct ExodusVariable {
    /// Variable name.
    pub name: String,
    /// Entity type: "node" or "element".
    pub entity_type: String,
    /// Variable values, one per entity.
    pub values: Vec<f64>,
}
/// Represents a boundary side set with element/side ID pairs and an optional
/// distance field (e.g. wall-normal distance).
pub struct ExodusSideSetData {
    /// Side-set identifier.
    pub id: u32,
    /// Human-readable name.
    pub name: String,
    /// Element IDs (1-based) in the set.
    pub element_ids: Vec<usize>,
    /// Local side IDs within each element.
    pub side_ids: Vec<usize>,
    /// Optional per-side scalar distribution factor.
    pub dist_factors: Vec<f64>,
}
/// Higher-level writing helpers for Exodus-II results.
pub struct ExodusWriter;
impl ExodusWriter {
    /// Append a time-step nodal variable record to a text buffer.
    ///
    /// The record format is:
    /// ```text
    /// NODAL_VAR `name` STEP `step` TIME `time`
    /// `v0` `v1` ... `vN`
    /// END_NODAL_VAR
    /// ```
    ///
    /// Returns a formatted `String`.
    pub fn write_nodal_variable_step(
        var_name: &str,
        step: usize,
        time: f64,
        values: &[f64],
    ) -> String {
        let mut s = String::new();
        s.push_str(&format!(
            "NODAL_VAR {} STEP {} TIME {:.8e}\n",
            var_name, step, time
        ));
        let vals: Vec<String> = values.iter().map(|v| format!("{:.8e}", v)).collect();
        s.push_str(&vals.join(" "));
        s.push('\n');
        s.push_str("END_NODAL_VAR\n");
        s
    }
    /// Serialize a global simulation variable to a one-line text record.
    ///
    /// Format:
    /// ```text
    /// GLOBAL_VAR `name` TIME `time` VALUE `value`
    /// ```
    pub fn write_global_variable(var: &ExodusGlobalVariable) -> String {
        format!(
            "GLOBAL_VAR {} TIME {:.8e} VALUE {:.8e}\n",
            var.name, var.time, var.value
        )
    }
    /// Serialize an `ExodusSideSetData` to a text block.
    ///
    /// Format:
    /// ```text
    /// SIDE_SET_DATA `id` `name` `count`
    /// `eid0` `sid0`
    /// `eid1` `sid1`
    /// ...
    /// [DIST_FACTORS
    /// `f0` `f1` ...]
    /// END_SIDE_SET_DATA
    /// ```
    pub fn write_side_set(data: &ExodusSideSetData) -> String {
        let mut s = String::new();
        let count = data.element_ids.len();
        s.push_str(&format!(
            "SIDE_SET_DATA {} {} {}\n",
            data.id, data.name, count
        ));
        for i in 0..count {
            s.push_str(&format!("{} {}\n", data.element_ids[i], data.side_ids[i]));
        }
        if !data.dist_factors.is_empty() {
            s.push_str("DIST_FACTORS\n");
            let fs: Vec<String> = data
                .dist_factors
                .iter()
                .map(|f| format!("{:.6e}", f))
                .collect();
            s.push_str(&fs.join(" "));
            s.push('\n');
        }
        s.push_str("END_SIDE_SET_DATA\n");
        s
    }
}
