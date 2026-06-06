//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::io::Write;

/// Dimension descriptor.
#[derive(Debug, Clone)]
pub struct NetCdfDimension {
    /// Dimension name.
    pub name: String,
    /// Dimension size. `None` indicates an unlimited dimension.
    pub size: Option<usize>,
}
impl NetCdfDimension {
    /// Check if this is an unlimited dimension.
    pub fn is_unlimited(&self) -> bool {
        self.size.is_none()
    }
    /// Get the effective size (0 for unlimited with no data).
    pub fn effective_size(&self) -> usize {
        self.size.unwrap_or(0)
    }
}
/// A single variable inside a [`NetCdfFile`].
#[derive(Debug, Clone)]
pub struct NetCdfVariable {
    /// Variable name.
    pub name: String,
    /// Ordered dimension names this variable depends on.
    pub dims: Vec<String>,
    /// Flat data values (row-major).
    pub data: Vec<f64>,
    /// Per-variable attributes.
    pub attributes: Vec<VariableAttribute>,
}
impl NetCdfVariable {
    /// Create a new variable with the given name, dimensions, and data.
    pub fn new(name: &str, dims: Vec<String>, data: Vec<f64>) -> Self {
        NetCdfVariable {
            name: name.to_string(),
            dims,
            data,
            attributes: Vec::new(),
        }
    }
    /// Add an attribute to this variable.
    pub fn add_attribute(&mut self, key: &str, value: &str) {
        self.attributes.push(VariableAttribute {
            key: key.to_string(),
            value: value.to_string(),
        });
    }
    /// Get the value of an attribute by key.
    pub fn get_attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|a| a.key == key)
            .map(|a| a.value.as_str())
    }
    /// Number of data elements.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Check if variable has no data.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
/// Extended statistics for a [`NetcdfVariable`].
#[derive(Debug, Clone)]
pub struct NetcdfVariableStats {
    /// Variable name.
    pub name: String,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// Standard deviation (population).
    pub std_dev: f64,
    /// Number of elements.
    pub count: usize,
}
impl NetcdfVariableStats {
    /// Range = max − min.
    pub fn range(&self) -> f64 {
        self.max - self.min
    }
    /// Coefficient of variation (std/mean), or 0 if mean == 0.
    pub fn cv(&self) -> f64 {
        if self.mean.abs() < f64::EPSILON {
            0.0
        } else {
            self.std_dev / self.mean
        }
    }
}
/// A NetCDF-4 (HDF5-based) in-memory dataset with group support.
#[derive(Debug, Clone)]
pub struct Nc4File {
    /// Root group (equivalent to "/" in HDF5).
    pub root: Nc4Group,
    /// Global (file-level) attributes.
    pub global_attributes: Vec<(String, String)>,
    /// All dimensions (name → size, None = unlimited).
    pub dimensions: Vec<(String, Option<usize>)>,
    /// Current size of the unlimited dimension (number of records written).
    pub unlimited_size: usize,
}
impl Nc4File {
    /// Create a new empty NetCDF-4 file.
    pub fn new() -> Self {
        Self {
            root: Nc4Group::new("/"),
            global_attributes: Vec::new(),
            dimensions: Vec::new(),
            unlimited_size: 0,
        }
    }
    /// Add a fixed dimension.
    pub fn add_dimension(&mut self, name: &str, size: usize) {
        self.dimensions.push((name.to_string(), Some(size)));
    }
    /// Add an unlimited dimension.
    pub fn add_unlimited_dimension(&mut self, name: &str) {
        self.dimensions.push((name.to_string(), None));
    }
    /// Extend the unlimited dimension by one record.
    pub fn extend_unlimited(&mut self) {
        self.unlimited_size += 1;
        for entry in self.dimensions.iter_mut() {
            if entry.1.is_none() || entry.1 == Some(self.unlimited_size - 1) {
                entry.1 = Some(self.unlimited_size);
                break;
            }
        }
    }
    /// Add a global attribute.
    pub fn add_global_attribute(&mut self, key: &str, value: &str) {
        self.global_attributes
            .push((key.to_string(), value.to_string()));
    }
    /// Get a global attribute.
    pub fn get_global_attribute(&self, key: &str) -> Option<&str> {
        self.global_attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
    /// Get a dimension size.
    pub fn get_dimension_size(&self, name: &str) -> Option<usize> {
        self.dimensions
            .iter()
            .find(|(n, _)| n == name)
            .and_then(|(_, s)| *s)
    }
    /// Whether a dimension is unlimited.
    pub fn is_unlimited(&self, name: &str) -> bool {
        self.dimensions.iter().any(|(n, s)| {
            n == name && s.is_none()
                || (n == name && *s == Some(self.unlimited_size) && self.unlimited_size > 0)
        })
    }
    /// Add a variable to the root group.
    pub fn add_variable(&mut self, var: Nc4Variable) {
        self.root.add_variable(var);
    }
    /// Get a variable from the root group.
    pub fn get_variable(&self, name: &str) -> Option<&Nc4Variable> {
        self.root.get_variable(name)
    }
    /// Add a named sub-group to the root.
    pub fn add_group(&mut self, group: Nc4Group) {
        self.root.add_subgroup(group);
    }
    /// Get a sub-group by name from the root.
    pub fn get_group(&self, name: &str) -> Option<&Nc4Group> {
        self.root.subgroups.iter().find(|g| g.name == name)
    }
    /// All variable count (all groups).
    pub fn total_variable_count(&self) -> usize {
        self.root.total_variable_count()
    }
    /// Serialize to a simple CDL-based binary envelope.
    ///
    /// Layout: `OXNC4` magic (5 bytes) + payload length u32 LE + CDL-like text payload.
    ///
    /// The CDL payload stores global attributes, dimensions, unlimited size,
    /// and root variable names as a simple key=value text block.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut lines = Vec::new();
        lines.push(format!("unlimited_size:{}", self.unlimited_size));
        for (k, v) in &self.global_attributes {
            lines.push(format!("attr:{}={}", k, v));
        }
        for (name, size) in &self.dimensions {
            match size {
                Some(s) => lines.push(format!("dim:{}={}", name, s)),
                None => lines.push(format!("dim:{}=UNLIMITED", name)),
            }
        }
        for var in &self.root.variables {
            lines.push(format!(
                "var:{}:{}:{}",
                var.name,
                var.dims.join(","),
                var.data_type.type_name()
            ));
        }
        let payload = lines.join("\n");
        let payload_bytes = payload.as_bytes();
        let mut buf = Vec::with_capacity(9 + payload_bytes.len());
        buf.extend_from_slice(b"OXNC4");
        buf.extend_from_slice(&(payload_bytes.len() as u32).to_le_bytes());
        buf.extend_from_slice(payload_bytes);
        buf
    }
}
impl Default for Nc4File {
    fn default() -> Self {
        Self::new()
    }
}
/// NetCDF4 variable type tags (mirrors HDF5 / NetCDF-4 data types).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Nc4DataType {
    /// 64-bit IEEE floating point.
    Float64,
    /// 32-bit IEEE floating point.
    Float32,
    /// 32-bit signed integer.
    Int32,
    /// 8-bit unsigned integer.
    UInt8,
    /// Variable-length string (NC_STRING).
    String,
}
impl Nc4DataType {
    /// Return the typical byte width (variable-length types return 0).
    pub fn byte_width(&self) -> usize {
        match self {
            Nc4DataType::Float64 => 8,
            Nc4DataType::Float32 => 4,
            Nc4DataType::Int32 => 4,
            Nc4DataType::UInt8 => 1,
            Nc4DataType::String => 0,
        }
    }
    /// NetCDF-4 type name string.
    pub fn type_name(&self) -> &'static str {
        match self {
            Nc4DataType::Float64 => "double",
            Nc4DataType::Float32 => "float",
            Nc4DataType::Int32 => "int",
            Nc4DataType::UInt8 => "ubyte",
            Nc4DataType::String => "string",
        }
    }
}
/// Fluent builder for writing NetCDF-convention trajectories.
#[derive(Debug, Clone, Default)]
pub struct NetcdfTrajectoryBuilder {
    /// Title / description of the trajectory.
    pub title: String,
    /// Application name that produced this trajectory.
    pub application: String,
    /// Frames in insertion order.
    pub frames: Vec<TrajectoryFrame>,
    /// Whether positions are stored in Å (true) or nm (false).
    pub use_angstroms: bool,
}
impl NetcdfTrajectoryBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        NetcdfTrajectoryBuilder {
            title: String::new(),
            application: "OxiPhysics".to_string(),
            frames: Vec::new(),
            use_angstroms: true,
        }
    }
    /// Set the trajectory title.
    pub fn with_title(mut self, title: &str) -> Self {
        self.title = title.to_string();
        self
    }
    /// Set the application name.
    pub fn with_application(mut self, app: &str) -> Self {
        self.application = app.to_string();
        self
    }
    /// Use Å as the position unit.
    pub fn in_angstroms(mut self) -> Self {
        self.use_angstroms = true;
        self
    }
    /// Use nm as the position unit.
    pub fn in_nanometres(mut self) -> Self {
        self.use_angstroms = false;
        self
    }
    /// Append a frame to the trajectory.
    pub fn add_frame(&mut self, frame: TrajectoryFrame) {
        self.frames.push(frame);
    }
    /// Number of frames stored.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }
    /// Number of atoms (taken from the first frame; 0 if empty).
    pub fn n_atoms(&self) -> usize {
        self.frames.first().map(|f| f.n_atoms()).unwrap_or(0)
    }
    /// Compute the RMSD trajectory (each frame vs. frame 0).
    pub fn rmsd_series(&self) -> Vec<f64> {
        if self.frames.is_empty() {
            return vec![];
        }
        let ref_frame = &self.frames[0];
        self.frames.iter().map(|f| f.rmsd_from(ref_frame)).collect()
    }
    /// Write a minimal CDL representation to a writer.
    pub fn write_cdl<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let n_atoms = self.n_atoms();
        let n_frames = self.frame_count();
        writeln!(writer, "netcdf trajectory {{")?;
        writeln!(writer, "dimensions:")?;
        writeln!(writer, "\tatom = {n_atoms} ;")?;
        writeln!(writer, "\tframe = UNLIMITED ; // currently {n_frames}")?;
        writeln!(writer, "\tspatial = 3 ;")?;
        writeln!(writer, "\tlabel = 5 ;")?;
        writeln!(writer, "variables:")?;
        let unit = if self.use_angstroms {
            "angstrom"
        } else {
            "nanometer"
        };
        writeln!(writer, "\tdouble coordinates(frame, atom, spatial) ;")?;
        writeln!(writer, "\t\tcoordinates:units = \"{unit}\" ;")?;
        writeln!(writer, "\tdouble time(frame) ;")?;
        writeln!(writer, "\t\ttime:units = \"picosecond\" ;")?;
        writeln!(writer, "\tdouble cell_lengths(frame, spatial) ;")?;
        writeln!(writer, "\t\tcell_lengths:units = \"{unit}\" ;")?;
        writeln!(writer, "// global attributes:")?;
        writeln!(writer, "\t:title = \"{}\" ;", self.title)?;
        writeln!(writer, "\t:application = \"{}\" ;", self.application)?;
        writeln!(writer, "\t:Conventions = \"AMBER\" ;")?;
        writeln!(writer, "\t:ConventionVersion = \"1.0\" ;")?;
        writeln!(writer, "}}")?;
        Ok(())
    }
    /// Extract the time series (ps) of all frames.
    pub fn time_series(&self) -> Vec<f64> {
        self.frames.iter().map(|f| f.time_ps).collect()
    }
    /// Compute per-frame centre-of-mass trajectories.
    pub fn com_trajectory(&self) -> Vec<[f64; 3]> {
        self.frames.iter().map(|f| f.centre_of_mass()).collect()
    }
    /// Get frame at a specific index (panics if out of range).
    pub fn frame(&self, idx: usize) -> &TrajectoryFrame {
        &self.frames[idx]
    }
}
/// An in-memory representation of a NetCDF-like dataset.
#[derive(Debug, Clone)]
pub struct NetCdfFile {
    /// Named dimensions and their sizes.
    pub dimensions: Vec<(String, usize)>,
    /// Variables stored in the file.
    pub variables: Vec<NetCdfVariable>,
    /// Global attributes as key-value string pairs.
    pub attributes: Vec<(String, String)>,
    /// Unlimited dimension name, if any.
    pub unlimited_dim: Option<String>,
}
impl NetCdfFile {
    /// Create an empty NetCDF file.
    pub fn new() -> Self {
        NetCdfFile {
            dimensions: Vec::new(),
            variables: Vec::new(),
            attributes: Vec::new(),
            unlimited_dim: None,
        }
    }
    /// Add a dimension.
    pub fn add_dimension(&mut self, name: &str, size: usize) {
        self.dimensions.push((name.to_string(), size));
    }
    /// Add an unlimited dimension.
    pub fn add_unlimited_dimension(&mut self, name: &str, current_size: usize) {
        self.dimensions.push((name.to_string(), current_size));
        self.unlimited_dim = Some(name.to_string());
    }
    /// Add a variable.
    pub fn add_variable(&mut self, var: NetCdfVariable) {
        self.variables.push(var);
    }
    /// Add a global attribute.
    pub fn add_attribute(&mut self, key: &str, value: &str) {
        self.attributes.push((key.to_string(), value.to_string()));
    }
    /// Get a variable by name.
    pub fn get_variable(&self, name: &str) -> Option<&NetCdfVariable> {
        self.variables.iter().find(|v| v.name == name)
    }
    /// Get a mutable variable by name.
    pub fn get_variable_mut(&mut self, name: &str) -> Option<&mut NetCdfVariable> {
        self.variables.iter_mut().find(|v| v.name == name)
    }
    /// Get a global attribute value by key.
    pub fn get_attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
    /// Get the size of a dimension by name.
    pub fn get_dimension_size(&self, name: &str) -> Option<usize> {
        self.dimensions
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, s)| *s)
    }
    /// Check if a dimension is unlimited.
    pub fn is_unlimited_dimension(&self, name: &str) -> bool {
        self.unlimited_dim.as_deref() == Some(name)
    }
    /// List all variable names.
    pub fn variable_names(&self) -> Vec<&str> {
        self.variables.iter().map(|v| v.name.as_str()).collect()
    }
    /// List all dimension names.
    pub fn dimension_names(&self) -> Vec<&str> {
        self.dimensions.iter().map(|(n, _)| n.as_str()).collect()
    }
}
impl Default for NetCdfFile {
    fn default() -> Self {
        Self::new()
    }
}
/// In-memory NetCDF-like file using HashMap for O(1) dimension lookups.
#[derive(Debug, Clone)]
pub struct NetcdfFile {
    /// Dimension name → size.
    pub dimensions: std::collections::HashMap<String, usize>,
    /// Variables.
    pub variables: Vec<NetcdfVariable>,
    /// Global attributes as (key, value) pairs.
    pub global_attrs: Vec<(String, String)>,
}
impl NetcdfFile {
    /// Create a new empty `NetcdfFile`.
    pub fn new() -> Self {
        Self {
            dimensions: std::collections::HashMap::new(),
            variables: Vec::new(),
            global_attrs: Vec::new(),
        }
    }
    /// Add a named dimension with size `size`.
    pub fn add_dimension(&mut self, name: &str, size: usize) {
        self.dimensions.insert(name.to_string(), size);
    }
    /// Add a variable.
    pub fn add_variable(&mut self, var: NetcdfVariable) {
        self.variables.push(var);
    }
    /// Add a global attribute.
    pub fn add_global_attr(&mut self, key: &str, value: &str) {
        self.global_attrs.push((key.to_string(), value.to_string()));
    }
    /// Serialize to a CDL text string.
    pub fn write_cdl(&self) -> String {
        let mut out = String::new();
        out.push_str("netcdf data {\n");
        out.push_str("dimensions:\n");
        let mut dims: Vec<(&String, &usize)> = self.dimensions.iter().collect();
        dims.sort_by_key(|(k, _)| k.as_str());
        for (name, size) in &dims {
            out.push_str(&format!("\t{} = {} ;\n", name, size));
        }
        out.push_str("variables:\n");
        for var in &self.variables {
            let dims_str = var.dimensions.join(", ");
            out.push_str(&format!("\tdouble {}({}) ;\n", var.name, dims_str));
            if !var.units.is_empty() {
                out.push_str(&format!("\t\t{}:units = \"{}\" ;\n", var.name, var.units));
            }
            if !var.long_name.is_empty() {
                out.push_str(&format!(
                    "\t\t{}:long_name = \"{}\" ;\n",
                    var.name, var.long_name
                ));
            }
        }
        if !self.global_attrs.is_empty() {
            out.push_str("// global attributes:\n");
            for (key, value) in &self.global_attrs {
                out.push_str(&format!("\t\t:{} = \"{}\" ;\n", key, value));
            }
        }
        out.push_str("data:\n");
        for var in &self.variables {
            let vals: Vec<String> = var.data.iter().map(|v| format!("{}", v)).collect();
            out.push_str(&format!("\t{} = {} ;\n", var.name, vals.join(", ")));
        }
        out.push_str("}\n");
        out
    }
    /// Write a simplified MD trajectory in NetCDF format (text CDL) to `path`.
    ///
    /// `times` is a slice of time values; `positions` is a slice of frames where
    /// each frame is a `Vec<[f64;3]>` of atom positions.
    pub fn trajectory_write(
        path: &str,
        times: &[f64],
        positions: &[Vec<[f64; 3]>],
    ) -> std::io::Result<()> {
        use std::io::Write;
        let n_frames = times.len();
        let n_atoms = positions.first().map(|f| f.len()).unwrap_or(0);
        let mut f = std::fs::File::create(path)?;
        writeln!(f, "netcdf trajectory {{")?;
        writeln!(f, "dimensions:")?;
        writeln!(f, "\tframe = {} ;", n_frames)?;
        writeln!(f, "\tatom = {} ;", n_atoms)?;
        writeln!(f, "\tspatial = 3 ;")?;
        writeln!(f, "variables:")?;
        writeln!(f, "\tdouble time(frame) ;")?;
        writeln!(f, "\t\ttime:units = \"ps\" ;")?;
        writeln!(f, "\tdouble coordinates(frame, atom, spatial) ;")?;
        writeln!(f, "\t\tcoordinates:units = \"angstrom\" ;")?;
        writeln!(f, "data:")?;
        let time_vals: Vec<String> = times.iter().map(|t| format!("{}", t)).collect();
        writeln!(f, "\ttime = {} ;", time_vals.join(", "))?;
        write!(f, "\tcoordinates = ")?;
        let mut first = true;
        for frame in positions {
            for pos in frame {
                for &coord in pos {
                    if !first {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", coord)?;
                    first = false;
                }
            }
        }
        writeln!(f, " ;")?;
        writeln!(f, "}}")?;
        Ok(())
    }
}
/// A builder-style writer for the simplified NetCDF-like text (CDL) format.
///
/// Supports dimension management, coordinate variables, typed variables,
/// global attributes, and trajectory serialisation.
#[derive(Debug, Clone)]
pub struct NetcdfWriter {
    pub(super) name: String,
    pub(super) dimensions: Vec<(String, usize)>,
    pub(super) unlimited_dim: Option<String>,
    pub(super) variables: Vec<NetcdfWriterVariable>,
    pub(super) global_attrs: Vec<(String, String)>,
}
impl NetcdfWriter {
    /// Create a new writer with a dataset name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            dimensions: Vec::new(),
            unlimited_dim: None,
            variables: Vec::new(),
            global_attrs: Vec::new(),
        }
    }
    /// Add a fixed-size dimension.
    pub fn add_dimension(&mut self, name: &str, size: usize) -> &mut Self {
        self.dimensions.push((name.to_string(), size));
        self
    }
    /// Add an unlimited (record) dimension.
    pub fn add_unlimited_dimension(&mut self, name: &str, current_size: usize) -> &mut Self {
        self.dimensions.push((name.to_string(), current_size));
        self.unlimited_dim = Some(name.to_string());
        self
    }
    /// Add a global attribute.
    pub fn add_global_attribute(&mut self, key: &str, value: &str) -> &mut Self {
        self.global_attrs.push((key.to_string(), value.to_string()));
        self
    }
    /// Add a variable with given dimension names and flat data.
    pub fn add_variable(&mut self, name: &str, dims: &[&str], data: Vec<f64>) -> &mut Self {
        self.variables.push(NetcdfWriterVariable {
            name: name.to_string(),
            dims: dims.iter().map(|s| s.to_string()).collect(),
            data,
            attrs: Vec::new(),
        });
        self
    }
    /// Add an attribute to the most recently added variable.
    pub fn add_variable_attribute(&mut self, key: &str, value: &str) -> &mut Self {
        if let Some(v) = self.variables.last_mut() {
            v.attrs.push((key.to_string(), value.to_string()));
        }
        self
    }
    /// Add a coordinate variable (1D variable with the same name as its dimension).
    pub fn add_coordinate(&mut self, dim_name: &str, data: Vec<f64>, units: &str) -> &mut Self {
        self.add_variable(dim_name, &[dim_name], data);
        self.add_variable_attribute("units", units);
        self.add_variable_attribute("axis", dim_name);
        self
    }
    /// Number of dimensions.
    pub fn n_dimensions(&self) -> usize {
        self.dimensions.len()
    }
    /// Number of variables.
    pub fn n_variables(&self) -> usize {
        self.variables.len()
    }
    /// Serialise to a CDL text string.
    pub fn to_cdl(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("netcdf {} {{\n", self.name));
        out.push_str("dimensions:\n");
        for (name, size) in &self.dimensions {
            if self.unlimited_dim.as_deref() == Some(name.as_str()) {
                out.push_str(&format!(
                    "\t{} = UNLIMITED ; // ({} currently)\n",
                    name, size
                ));
            } else {
                out.push_str(&format!("\t{} = {} ;\n", name, size));
            }
        }
        out.push_str("variables:\n");
        for var in &self.variables {
            let dims = var.dims.join(", ");
            out.push_str(&format!("\tdouble {}({}) ;\n", var.name, dims));
            for (k, v) in &var.attrs {
                out.push_str(&format!("\t\t{}:{} = \"{}\" ;\n", var.name, k, v));
            }
        }
        if !self.global_attrs.is_empty() {
            out.push_str("// global attributes:\n");
            for (k, v) in &self.global_attrs {
                out.push_str(&format!("\t:{} = \"{}\" ;\n", k, v));
            }
        }
        out.push_str("data:\n");
        for var in &self.variables {
            let vals: Vec<String> = var.data.iter().map(|v| format!("{}", v)).collect();
            out.push_str(&format!("\t{} = {} ;\n", var.name, vals.join(", ")));
        }
        out.push_str("}\n");
        out
    }
    /// Consume the writer and return the CDL string.
    pub fn finish(self) -> String {
        self.to_cdl()
    }
    /// Write the CDL to a file.
    pub fn write_to_file(&self, path: &str) -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        f.write_all(self.to_cdl().as_bytes())
    }
}
/// A single variable in the new simplified NetcdfFile.
#[derive(Debug, Clone)]
pub struct NetcdfVariable {
    /// Variable name.
    pub name: String,
    /// Dimension names for this variable.
    pub dimensions: Vec<String>,
    /// Flat data values.
    pub data: Vec<f64>,
    /// Unit string (e.g., `"m/s"`).
    pub units: String,
    /// Long descriptive name.
    pub long_name: String,
}
/// One snapshot of positions + optional velocities for an MD trajectory.
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Simulation time in ps.
    pub time_ps: f64,
    /// Positions: `n_atoms × 3`, row-major (Å or nm depending on convention).
    pub positions: Vec<[f64; 3]>,
    /// Optional velocities: `n_atoms × 3` in nm/ps.
    pub velocities: Option<Vec<[f64; 3]>>,
    /// Optional box vectors (orthorhombic: `[Lx, Ly, Lz]`).
    pub box_lengths: Option<[f64; 3]>,
}
impl TrajectoryFrame {
    /// Number of atoms in this frame.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }
    /// Compute the centre-of-mass position (equal masses assumed).
    pub fn centre_of_mass(&self) -> [f64; 3] {
        if self.positions.is_empty() {
            return [0.0; 3];
        }
        let mut s = [0.0_f64; 3];
        for p in &self.positions {
            s[0] += p[0];
            s[1] += p[1];
            s[2] += p[2];
        }
        let inv = 1.0 / self.positions.len() as f64;
        [s[0] * inv, s[1] * inv, s[2] * inv]
    }
    /// Compute the root-mean-square displacement from a reference frame.
    pub fn rmsd_from(&self, reference: &TrajectoryFrame) -> f64 {
        let n = self.positions.len().min(reference.positions.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n)
            .map(|i| {
                let dx = self.positions[i][0] - reference.positions[i][0];
                let dy = self.positions[i][1] - reference.positions[i][1];
                let dz = self.positions[i][2] - reference.positions[i][2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum / n as f64).sqrt()
    }
    /// Compute kinetic energy assuming all masses equal to `mass_amu`.
    /// Returns energy in kJ/mol (1 amu × nm²/ps² = 1 kJ/mol).
    pub fn kinetic_energy(&self, mass_amu: f64) -> f64 {
        let vels = match &self.velocities {
            Some(v) => v,
            None => return 0.0,
        };
        let sum: f64 = vels
            .iter()
            .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum();
        0.5 * mass_amu * sum
    }
    /// Translate all positions by `delta`.
    pub fn translate(&mut self, delta: [f64; 3]) {
        for p in &mut self.positions {
            p[0] += delta[0];
            p[1] += delta[1];
            p[2] += delta[2];
        }
    }
}
/// A simple CDL text reader.
///
/// Reads dimension, variable, and data sections from CDL output produced by
/// [`NetcdfWriter`] or [`NetcdfFile::write_cdl`].
#[derive(Debug, Clone)]
pub struct NetcdfReader {
    /// Parsed dimensions (name → size).
    pub dimensions: std::collections::HashMap<String, usize>,
    /// Parsed variables.
    pub variables: Vec<NetcdfVariable>,
    /// Global attributes.
    pub global_attrs: Vec<(String, String)>,
    /// Unlimited dimension name.
    pub unlimited_dim: Option<String>,
}
impl NetcdfReader {
    /// Parse a CDL string.
    pub fn from_cdl(cdl: &str) -> Result<Self, String> {
        let file = parse_ncf_text(cdl)?;
        let variables: Vec<NetcdfVariable> = file
            .variables
            .into_iter()
            .map(|v| {
                let units = v.get_attribute("units").unwrap_or("").to_string();
                let long_name = v.get_attribute("long_name").unwrap_or("").to_string();
                NetcdfVariable {
                    name: v.name,
                    dimensions: v.dims,
                    data: v.data,
                    units,
                    long_name,
                }
            })
            .collect();
        Ok(Self {
            dimensions: file.dimensions.into_iter().collect(),
            variables,
            global_attrs: file.attributes,
            unlimited_dim: file.unlimited_dim,
        })
    }
    /// Get dimension size by name.
    pub fn get_dimension(&self, name: &str) -> Option<usize> {
        self.dimensions.get(name).copied()
    }
    /// Get variable by name.
    pub fn get_variable(&self, name: &str) -> Option<&NetcdfVariable> {
        self.variables.iter().find(|v| v.name == name)
    }
    /// Get data slice for a variable.
    pub fn get_data(&self, var_name: &str) -> Option<&[f64]> {
        self.get_variable(var_name).map(|v| v.data.as_slice())
    }
    /// List all variable names.
    pub fn variable_names(&self) -> Vec<&str> {
        self.variables.iter().map(|v| v.name.as_str()).collect()
    }
    /// List all dimension names.
    pub fn dimension_names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.dimensions.keys().map(String::as_str).collect();
        names.sort();
        names
    }
}
/// A variable attribute (key-value string pair).
#[derive(Debug, Clone)]
pub struct VariableAttribute {
    /// Attribute key.
    pub key: String,
    /// Attribute value.
    pub value: String,
}
/// A named group within a NetCDF-4 file.
///
/// Groups allow hierarchical organization of variables, akin to HDF5 groups.
#[derive(Debug, Clone)]
pub struct Nc4Group {
    /// Group name.
    pub name: String,
    /// Variables in this group.
    pub variables: Vec<Nc4Variable>,
    /// Sub-groups.
    pub subgroups: Vec<Nc4Group>,
    /// Group-level attributes (inherited by sub-groups).
    pub attributes: Vec<(String, String)>,
}
impl Nc4Group {
    /// Create a new empty group.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            variables: Vec::new(),
            subgroups: Vec::new(),
            attributes: Vec::new(),
        }
    }
    /// Add a variable to this group.
    pub fn add_variable(&mut self, var: Nc4Variable) {
        self.variables.push(var);
    }
    /// Add a sub-group.
    pub fn add_subgroup(&mut self, group: Nc4Group) {
        self.subgroups.push(group);
    }
    /// Add an attribute.
    pub fn add_attribute(&mut self, key: &str, value: &str) {
        self.attributes.push((key.to_string(), value.to_string()));
    }
    /// Get a variable by name.
    pub fn get_variable(&self, name: &str) -> Option<&Nc4Variable> {
        self.variables.iter().find(|v| v.name == name)
    }
    /// Get attribute value.
    pub fn get_attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
    /// All variables in this group and all sub-groups (depth-first).
    pub fn all_variables(&self) -> Vec<&Nc4Variable> {
        let mut out: Vec<&Nc4Variable> = self.variables.iter().collect();
        for sg in &self.subgroups {
            out.extend(sg.all_variables());
        }
        out
    }
    /// Total variable count (recursive).
    pub fn total_variable_count(&self) -> usize {
        self.variables.len()
            + self
                .subgroups
                .iter()
                .map(|g| g.total_variable_count())
                .sum::<usize>()
    }
}
/// A dimension with name, size, and whether it is an unlimited record dimension.
#[derive(Debug, Clone, PartialEq)]
pub struct NetcdfDimSpec {
    /// Dimension name.
    pub name: String,
    /// Current size (may grow if unlimited).
    pub size: usize,
    /// Whether this is the unlimited (record) dimension.
    pub unlimited: bool,
}
impl NetcdfDimSpec {
    /// Create a fixed-size dimension.
    pub fn fixed(name: &str, size: usize) -> Self {
        NetcdfDimSpec {
            name: name.to_string(),
            size,
            unlimited: false,
        }
    }
    /// Create an unlimited record dimension with current size.
    pub fn unlimited(name: &str, current_size: usize) -> Self {
        NetcdfDimSpec {
            name: name.to_string(),
            size: current_size,
            unlimited: true,
        }
    }
    /// Produce the CDL declaration string.
    pub fn to_cdl(&self) -> String {
        if self.unlimited {
            format!("\t{} = UNLIMITED ; // currently {}", self.name, self.size)
        } else {
            format!("\t{} = {} ;", self.name, self.size)
        }
    }
}
/// A NetCDF-4 variable with typed storage.
#[derive(Debug, Clone)]
pub struct Nc4Variable {
    /// Variable name.
    pub name: String,
    /// Ordered dimension names.
    pub dims: Vec<String>,
    /// Floating-point data (used when data_type is Float64 or Float32).
    pub float_data: Vec<f64>,
    /// String data (used when data_type is String).
    pub string_data: Vec<String>,
    /// Integer data (used when data_type is Int32 or UInt8).
    pub int_data: Vec<i64>,
    /// Data type.
    pub data_type: Nc4DataType,
    /// Per-variable attributes.
    pub attributes: Vec<VariableAttribute>,
}
impl Nc4Variable {
    /// Create a Float64 variable.
    pub fn float64(name: &str, dims: Vec<String>, data: Vec<f64>) -> Self {
        Self {
            name: name.to_string(),
            dims,
            float_data: data,
            string_data: Vec::new(),
            int_data: Vec::new(),
            data_type: Nc4DataType::Float64,
            attributes: Vec::new(),
        }
    }
    /// Create a String variable.
    pub fn string_var(name: &str, dims: Vec<String>, data: Vec<String>) -> Self {
        Self {
            name: name.to_string(),
            dims,
            float_data: Vec::new(),
            string_data: data,
            int_data: Vec::new(),
            data_type: Nc4DataType::String,
            attributes: Vec::new(),
        }
    }
    /// Create an Int32 variable.
    pub fn int32(name: &str, dims: Vec<String>, data: Vec<i64>) -> Self {
        Self {
            name: name.to_string(),
            dims,
            float_data: Vec::new(),
            string_data: Vec::new(),
            int_data: data,
            data_type: Nc4DataType::Int32,
            attributes: Vec::new(),
        }
    }
    /// Add an attribute.
    pub fn add_attribute(&mut self, key: &str, value: &str) {
        self.attributes.push(VariableAttribute {
            key: key.to_string(),
            value: value.to_string(),
        });
    }
    /// Get attribute value by key.
    pub fn get_attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|a| a.key == key)
            .map(|a| a.value.as_str())
    }
    /// Number of elements (type-appropriate).
    pub fn len(&self) -> usize {
        match self.data_type {
            Nc4DataType::String => self.string_data.len(),
            Nc4DataType::Int32 | Nc4DataType::UInt8 => self.int_data.len(),
            _ => self.float_data.len(),
        }
    }
    /// Whether the variable has no data.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
/// Internal variable representation for [`NetcdfWriter`].
#[derive(Debug, Clone)]
pub(super) struct NetcdfWriterVariable {
    pub(super) name: String,
    pub(super) dims: Vec<String>,
    pub(super) data: Vec<f64>,
    pub(super) attrs: Vec<(String, String)>,
}
