//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{foam_header, parse_dict_tokens, strip_foam_comments, tokenise_foam};

/// A single boundary condition entry for a `FoamField`.
pub struct FoamBc {
    /// Patch name matching an entry in the boundary file.
    pub patch_name: String,
    /// Boundary condition type: "zeroGradient", "fixedValue", "noSlip", etc.
    pub bc_type: String,
    /// Optional value string used with "fixedValue" conditions.
    pub value: Option<String>,
}
/// Common boundary condition constructors.
impl FoamBc {
    /// Zero-gradient (Neumann) condition.
    pub fn zero_gradient(patch: &str) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "zeroGradient".to_string(),
            value: None,
        }
    }
    /// Fixed scalar value (Dirichlet) condition.
    pub fn fixed_scalar(patch: &str, val: f64) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "fixedValue".to_string(),
            value: Some(format!("uniform {}", val)),
        }
    }
    /// Fixed vector value (Dirichlet) condition.
    pub fn fixed_vector(patch: &str, val: [f64; 3]) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "fixedValue".to_string(),
            value: Some(format!("uniform ({} {} {})", val[0], val[1], val[2])),
        }
    }
    /// No-slip wall condition.
    pub fn no_slip(patch: &str) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "noSlip".to_string(),
            value: None,
        }
    }
    /// Symmetry condition.
    pub fn symmetry(patch: &str) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "symmetry".to_string(),
            value: None,
        }
    }
    /// Inlet-outlet condition with a given inlet value.
    pub fn inlet_outlet_scalar(patch: &str, inlet_val: f64) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "inletOutlet".to_string(),
            value: Some(format!("uniform {}", inlet_val)),
        }
    }
    /// Pressure inlet-outlet condition.
    pub fn pressure_inlet_outlet(patch: &str, val: f64) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "totalPressure".to_string(),
            value: Some(format!("uniform {}", val)),
        }
    }
    /// Empty condition (for 2D cases).
    pub fn empty(patch: &str) -> Self {
        FoamBc {
            patch_name: patch.to_string(),
            bc_type: "empty".to_string(),
            value: None,
        }
    }
}
/// A polyhedral mesh for OpenFOAM output.
pub struct FoamMesh {
    /// Vertex coordinates.
    pub points: Vec<[f64; 3]>,
    /// Face vertex lists (each face is an arbitrary polygon).
    pub faces: Vec<Vec<usize>>,
    /// Owner cell index for each face.
    pub owner: Vec<usize>,
    /// Neighbour cell index for each face (-1 for boundary faces).
    pub neighbour: Vec<i64>,
    /// Total number of cells.
    pub n_cells: usize,
    /// Boundary patch definitions.
    pub boundary_patches: Vec<FoamPatch>,
}
impl FoamMesh {
    /// Build a simple structured hex mesh: `nx × ny × nz` cells filling
    /// a box `[0, lx] × [0, ly] × [0, lz]`.
    ///
    /// Face ordering follows the OpenFOAM convention (lower-index owner,
    /// internal faces first, then boundary faces grouped by patch).
    pub fn box_mesh(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> Self {
        let mut points = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    points.push([
                        lx * i as f64 / nx as f64,
                        ly * j as f64 / ny as f64,
                        lz * k as f64 / nz as f64,
                    ]);
                }
            }
        }
        let pid =
            |i: usize, j: usize, k: usize| -> usize { k * (ny + 1) * (nx + 1) + j * (nx + 1) + i };
        let cid = |i: usize, j: usize, k: usize| -> usize { k * ny * nx + j * nx + i };
        let n_cells = nx * ny * nz;
        let mut faces: Vec<Vec<usize>> = Vec::new();
        let mut owner: Vec<usize> = Vec::new();
        let mut neighbour: Vec<i64> = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 1..nx {
                    faces.push(vec![
                        pid(i, j, k),
                        pid(i, j + 1, k),
                        pid(i, j + 1, k + 1),
                        pid(i, j, k + 1),
                    ]);
                    owner.push(cid(i - 1, j, k));
                    neighbour.push(cid(i, j, k) as i64);
                }
            }
        }
        for k in 0..nz {
            for j in 1..ny {
                for i in 0..nx {
                    faces.push(vec![
                        pid(i, j, k),
                        pid(i + 1, j, k),
                        pid(i + 1, j, k + 1),
                        pid(i, j, k + 1),
                    ]);
                    owner.push(cid(i, j - 1, k));
                    neighbour.push(cid(i, j, k) as i64);
                }
            }
        }
        for k in 1..nz {
            for j in 0..ny {
                for i in 0..nx {
                    faces.push(vec![
                        pid(i, j, k),
                        pid(i + 1, j, k),
                        pid(i + 1, j + 1, k),
                        pid(i, j + 1, k),
                    ]);
                    owner.push(cid(i, j, k - 1));
                    neighbour.push(cid(i, j, k) as i64);
                }
            }
        }
        let xmin_start = faces.len();
        for k in 0..nz {
            for j in 0..ny {
                faces.push(vec![
                    pid(0, j, k),
                    pid(0, j, k + 1),
                    pid(0, j + 1, k + 1),
                    pid(0, j + 1, k),
                ]);
                owner.push(cid(0, j, k));
                neighbour.push(-1);
            }
        }
        let xmin_n = faces.len() - xmin_start;
        let xmax_start = faces.len();
        for k in 0..nz {
            for j in 0..ny {
                faces.push(vec![
                    pid(nx, j, k),
                    pid(nx, j + 1, k),
                    pid(nx, j + 1, k + 1),
                    pid(nx, j, k + 1),
                ]);
                owner.push(cid(nx - 1, j, k));
                neighbour.push(-1);
            }
        }
        let xmax_n = faces.len() - xmax_start;
        let ymin_start = faces.len();
        for k in 0..nz {
            for i in 0..nx {
                faces.push(vec![
                    pid(i, 0, k),
                    pid(i + 1, 0, k),
                    pid(i + 1, 0, k + 1),
                    pid(i, 0, k + 1),
                ]);
                owner.push(cid(i, 0, k));
                neighbour.push(-1);
            }
        }
        let ymin_n = faces.len() - ymin_start;
        let ymax_start = faces.len();
        for k in 0..nz {
            for i in 0..nx {
                faces.push(vec![
                    pid(i, ny, k),
                    pid(i, ny, k + 1),
                    pid(i + 1, ny, k + 1),
                    pid(i + 1, ny, k),
                ]);
                owner.push(cid(i, ny - 1, k));
                neighbour.push(-1);
            }
        }
        let ymax_n = faces.len() - ymax_start;
        let zmin_start = faces.len();
        for j in 0..ny {
            for i in 0..nx {
                faces.push(vec![
                    pid(i, j, 0),
                    pid(i, j + 1, 0),
                    pid(i + 1, j + 1, 0),
                    pid(i + 1, j, 0),
                ]);
                owner.push(cid(i, j, 0));
                neighbour.push(-1);
            }
        }
        let zmin_n = faces.len() - zmin_start;
        let zmax_start = faces.len();
        for j in 0..ny {
            for i in 0..nx {
                faces.push(vec![
                    pid(i, j, nz),
                    pid(i + 1, j, nz),
                    pid(i + 1, j + 1, nz),
                    pid(i, j + 1, nz),
                ]);
                owner.push(cid(i, j, nz - 1));
                neighbour.push(-1);
            }
        }
        let zmax_n = faces.len() - zmax_start;
        let boundary_patches = vec![
            FoamPatch {
                name: "xmin".into(),
                patch_type: "patch".into(),
                start_face: xmin_start,
                n_faces: xmin_n,
            },
            FoamPatch {
                name: "xmax".into(),
                patch_type: "patch".into(),
                start_face: xmax_start,
                n_faces: xmax_n,
            },
            FoamPatch {
                name: "ymin".into(),
                patch_type: "patch".into(),
                start_face: ymin_start,
                n_faces: ymin_n,
            },
            FoamPatch {
                name: "ymax".into(),
                patch_type: "patch".into(),
                start_face: ymax_start,
                n_faces: ymax_n,
            },
            FoamPatch {
                name: "zmin".into(),
                patch_type: "patch".into(),
                start_face: zmin_start,
                n_faces: zmin_n,
            },
            FoamPatch {
                name: "zmax".into(),
                patch_type: "patch".into(),
                start_face: zmax_start,
                n_faces: zmax_n,
            },
        ];
        FoamMesh {
            points,
            faces,
            owner,
            neighbour,
            n_cells,
            boundary_patches,
        }
    }
    /// Write the `points` file content as a `String`.
    pub fn write_points(&self) -> String {
        let mut s = foam_header("vectorField", "points");
        s.push('\n');
        s.push_str(&format!("{}\n(\n", self.points.len()));
        for p in &self.points {
            s.push_str(&format!("({} {} {})\n", p[0], p[1], p[2]));
        }
        s.push_str(")\n");
        s
    }
    /// Write the `faces` file content as a `String`.
    pub fn write_faces(&self) -> String {
        let mut s = foam_header("faceList", "faces");
        s.push('\n');
        s.push_str(&format!("{}\n(\n", self.faces.len()));
        for face in &self.faces {
            let verts: Vec<String> = face.iter().map(|v| v.to_string()).collect();
            s.push_str(&format!("{}({})\n", face.len(), verts.join(" ")));
        }
        s.push_str(")\n");
        s
    }
    /// Write the `owner` file content as a `String`.
    pub fn write_owner(&self) -> String {
        let mut s = foam_header("labelList", "owner");
        s.push('\n');
        s.push_str(&format!("{}\n(\n", self.owner.len()));
        for &o in &self.owner {
            s.push_str(&format!("{}\n", o));
        }
        s.push_str(")\n");
        s
    }
    /// Write the `neighbour` file content as a `String`.
    pub fn write_neighbour(&self) -> String {
        let mut s = foam_header("labelList", "neighbour");
        s.push('\n');
        s.push_str(&format!("{}\n(\n", self.neighbour.len()));
        for &n in &self.neighbour {
            s.push_str(&format!("{}\n", n));
        }
        s.push_str(")\n");
        s
    }
    /// Write the `boundary` file content as a `String`.
    pub fn write_boundary(&self) -> String {
        let mut s = foam_header("polyBoundaryMesh", "boundary");
        s.push('\n');
        s.push_str(&format!("{}\n(\n", self.boundary_patches.len()));
        for patch in &self.boundary_patches {
            s.push_str(
                &format!(
                    "    {}\n    {{\n        type        {};\n        nFaces      {};\n        startFace   {};\n    }}\n",
                    patch.name, patch.patch_type, patch.n_faces, patch.start_face
                ),
            );
        }
        s.push_str(")\n");
        s
    }
}
impl FoamMesh {
    /// Count internal faces (those with a non-negative neighbour).
    pub fn n_internal_faces(&self) -> usize {
        self.neighbour.iter().filter(|&&n| n >= 0).count()
    }
    /// Count boundary faces (those with neighbour == -1).
    pub fn n_boundary_faces(&self) -> usize {
        self.neighbour.iter().filter(|&&n| n < 0).count()
    }
    /// Total number of faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }
    /// Total number of points.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
    /// Compute the bounding box of the mesh: `([xmin, ymin, zmin], [xmax, ymax, zmax])`.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.points.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut min = self.points[0];
        let mut max = self.points[0];
        for p in &self.points {
            for i in 0..3 {
                if p[i] < min[i] {
                    min[i] = p[i];
                }
                if p[i] > max[i] {
                    max[i] = p[i];
                }
            }
        }
        (min, max)
    }
    /// Compute the centre of the mesh bounding box.
    pub fn centre(&self) -> [f64; 3] {
        let (min, max) = self.bounding_box();
        [
            (min[0] + max[0]) * 0.5,
            (min[1] + max[1]) * 0.5,
            (min[2] + max[2]) * 0.5,
        ]
    }
    /// Check basic mesh consistency: owner/neighbour length matches faces.
    pub fn check_topology(&self) -> bool {
        self.owner.len() == self.faces.len() && self.neighbour.len() == self.faces.len()
    }
    /// Find the patch with a given name.
    pub fn find_patch(&self, name: &str) -> Option<&FoamPatch> {
        self.boundary_patches.iter().find(|p| p.name == name)
    }
    /// Return all patch names.
    pub fn patch_names(&self) -> Vec<&str> {
        self.boundary_patches
            .iter()
            .map(|p| p.name.as_str())
            .collect()
    }
}
/// Represents a single time directory's metadata.
#[derive(Debug, Clone)]
pub struct FoamTimeDir {
    /// Time value (e.g. 0.0, 0.5, 1.0).
    pub time: f64,
    /// Directory name as a string (e.g. "0", "0.5", "1").
    pub dir_name: String,
    /// List of field names found at this time (e.g. \["U", "p"\]).
    pub fields: Vec<String>,
}
/// Writer for the OpenFOAM `system/controlDict` file.
pub struct ControlDict {
    /// Solver application name (e.g. "icoFoam", "simpleFoam").
    pub application: String,
    /// Simulation start time.
    pub start_time: f64,
    /// Simulation end time.
    pub end_time: f64,
    /// Time step size.
    pub delta_t: f64,
    /// How often to write results (in time units).
    pub write_interval: f64,
    /// Output format: "ascii" or "binary".
    pub write_format: String,
    /// Number of significant digits for ASCII output.
    pub write_precision: usize,
}
impl ControlDict {
    /// Create a `ControlDict` with sensible defaults.
    pub fn new(application: &str, end_time: f64, dt: f64) -> Self {
        ControlDict {
            application: application.to_string(),
            start_time: 0.0,
            end_time,
            delta_t: dt,
            write_interval: end_time / 10.0,
            write_format: "ascii".to_string(),
            write_precision: 6,
        }
    }
}
impl std::fmt::Display for ControlDict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = foam_header("dictionary", "controlDict");
        s.push('\n');
        s.push_str(&format!("application     {};\n\n", self.application));
        s.push_str("startFrom       startTime;\n\n");
        s.push_str(&format!("startTime       {};\n\n", self.start_time));
        s.push_str("stopAt          endTime;\n\n");
        s.push_str(&format!("endTime         {};\n\n", self.end_time));
        s.push_str(&format!("deltaT          {};\n\n", self.delta_t));
        s.push_str(&format!("writeFormat     {};\n\n", self.write_format));
        s.push_str(&format!("writePrecision  {};\n\n", self.write_precision));
        s.push_str("writeCompression off;\n\n");
        s.push_str("timeFormat      general;\n\n");
        s.push_str("timePrecision   6;\n\n");
        s.push_str("runTimeModifiable true;\n\n");
        s.push_str(&format!("writeInterval   {};\n", self.write_interval));
        write!(f, "{s}")
    }
}
/// Boundary patch definition for OpenFOAM output.
#[derive(Debug, Clone)]
pub struct FoamPatch {
    /// Name of the patch (e.g. "inlet", "outlet", "walls").
    pub name: String,
    /// Patch type: "wall", "patch", "symmetryPlane", etc.
    pub patch_type: String,
    /// Index of the first face belonging to this patch.
    pub start_face: usize,
    /// Number of faces in this patch.
    pub n_faces: usize,
}
/// A parsed OpenFOAM dictionary (key-value map preserving insertion order).
#[derive(Debug, Clone, PartialEq)]
pub struct FoamDict {
    /// Ordered key-value entries.
    pub entries: Vec<(String, FoamValue)>,
}
impl FoamDict {
    /// Create an empty dictionary.
    pub fn new() -> Self {
        FoamDict {
            entries: Vec::new(),
        }
    }
    /// Insert a key-value pair (appends; does not deduplicate).
    pub fn insert(&mut self, key: impl Into<String>, value: FoamValue) {
        self.entries.push((key.into(), value));
    }
    /// Look up a key and return the first matching value.
    pub fn get(&self, key: &str) -> Option<&FoamValue> {
        self.entries.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }
    /// Look up a key and return it as a scalar, if possible.
    pub fn get_scalar(&self, key: &str) -> Option<f64> {
        match self.get(key) {
            Some(FoamValue::Scalar(v)) => Some(*v),
            _ => None,
        }
    }
    /// Look up a key and return it as a word/string, if possible.
    pub fn get_word(&self, key: &str) -> Option<&str> {
        match self.get(key) {
            Some(FoamValue::Word(w)) => Some(w.as_str()),
            _ => None,
        }
    }
    /// Look up a key and return it as a sub-dictionary, if possible.
    pub fn get_dict(&self, key: &str) -> Option<&FoamDict> {
        match self.get(key) {
            Some(FoamValue::Dict(d)) => Some(d),
            _ => None,
        }
    }
    /// Look up a key and return it as a vector, if possible.
    pub fn get_vector(&self, key: &str) -> Option<[f64; 3]> {
        match self.get(key) {
            Some(FoamValue::Vector(v)) => Some(*v),
            _ => None,
        }
    }
    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the dictionary is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// All keys in insertion order.
    pub fn keys(&self) -> Vec<&str> {
        self.entries.iter().map(|(k, _)| k.as_str()).collect()
    }
    /// Serialise to OpenFOAM dictionary format string.
    pub fn to_foam_string(&self, indent: usize) -> String {
        let pad = "    ".repeat(indent);
        let mut s = String::new();
        for (key, value) in &self.entries {
            match value {
                FoamValue::Scalar(v) => {
                    s.push_str(&format!("{}{:<16}{};\n", pad, key, v));
                }
                FoamValue::Word(w) => {
                    s.push_str(&format!("{}{:<16}{};\n", pad, key, w));
                }
                FoamValue::Vector(v) => {
                    s.push_str(&format!(
                        "{}{:<16}({} {} {});\n",
                        pad, key, v[0], v[1], v[2]
                    ));
                }
                FoamValue::Dict(d) => {
                    s.push_str(&format!("{}{}\n{}{{\n", pad, key, pad));
                    s.push_str(&d.to_foam_string(indent + 1));
                    s.push_str(&format!("{}}}\n", pad));
                }
                FoamValue::List(items) => {
                    s.push_str(&format!("{}{}\n{}(\n", pad, key, pad));
                    for item in items {
                        match item {
                            FoamValue::Scalar(v) => {
                                s.push_str(&format!("{}    {}\n", pad, v));
                            }
                            FoamValue::Word(w) => {
                                s.push_str(&format!("{}    {}\n", pad, w));
                            }
                            FoamValue::Vector(v) => {
                                s.push_str(&format!("{}    ({} {} {})\n", pad, v[0], v[1], v[2]));
                            }
                            _ => {}
                        }
                    }
                    s.push_str(&format!("{});\n", pad));
                }
            }
        }
        s
    }
    /// Parse a simplified OpenFOAM dictionary from a string.
    ///
    /// Supports: `key value;`, `key { ... }`, `key ( ... );`, and
    /// `key (x y z);` vector syntax.  Comments (`//` and `/* */`) are stripped.
    pub fn parse(input: &str) -> Self {
        let cleaned = strip_foam_comments(input);
        let tokens = tokenise_foam(&cleaned);
        let (dict, _) = parse_dict_tokens(&tokens, 0);
        dict
    }
}
impl Default for FoamDict {
    fn default() -> Self {
        Self::new()
    }
}
/// Writer for the OpenFOAM `system/fvSchemes` file.
pub struct FvSchemes {
    /// Time derivative scheme.
    pub ddt_scheme: String,
    /// Gradient scheme.
    pub grad_scheme: String,
    /// Divergence schemes (key → scheme).
    pub div_schemes: Vec<(String, String)>,
    /// Laplacian schemes (key → scheme).
    pub laplacian_schemes: Vec<(String, String)>,
    /// Interpolation scheme.
    pub interpolation_scheme: String,
    /// Surface-normal gradient scheme.
    pub sn_grad_scheme: String,
}
impl FvSchemes {
    /// Create default second-order schemes.
    pub fn default_second_order() -> Self {
        FvSchemes {
            ddt_scheme: "Euler".to_string(),
            grad_scheme: "Gauss linear".to_string(),
            div_schemes: vec![
                ("default".to_string(), "none".to_string()),
                (
                    "div(phi,U)".to_string(),
                    "Gauss linearUpwind grad(U)".to_string(),
                ),
            ],
            laplacian_schemes: vec![("default".to_string(), "Gauss linear corrected".to_string())],
            interpolation_scheme: "linear".to_string(),
            sn_grad_scheme: "corrected".to_string(),
        }
    }
}
impl std::fmt::Display for FvSchemes {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = foam_header("dictionary", "fvSchemes");
        s.push('\n');
        s.push_str("ddtSchemes\n{\n");
        s.push_str(&format!("    default         {};\n", self.ddt_scheme));
        s.push_str("}\n\n");
        s.push_str("gradSchemes\n{\n");
        s.push_str(&format!("    default         {};\n", self.grad_scheme));
        s.push_str("}\n\n");
        s.push_str("divSchemes\n{\n");
        for (key, val) in &self.div_schemes {
            s.push_str(&format!("    {}    {};\n", key, val));
        }
        s.push_str("}\n\n");
        s.push_str("laplacianSchemes\n{\n");
        for (key, val) in &self.laplacian_schemes {
            s.push_str(&format!("    {}    {};\n", key, val));
        }
        s.push_str("}\n\n");
        s.push_str("interpolationSchemes\n{\n");
        s.push_str(&format!(
            "    default         {};\n",
            self.interpolation_scheme
        ));
        s.push_str("}\n\n");
        s.push_str("snGradSchemes\n{\n");
        s.push_str(&format!("    default         {};\n", self.sn_grad_scheme));
        s.push_str("}\n");
        write!(f, "{s}")
    }
}
/// A parsed OpenFOAM solver residual entry (from log file).
#[derive(Debug, Clone)]
pub struct FoamResidual {
    /// Time step (or iteration) number.
    pub time: f64,
    /// Field name (e.g., "p", "U").
    pub field: String,
    /// Initial residual.
    pub initial_residual: f64,
    /// Final residual.
    pub final_residual: f64,
    /// Number of iterations.
    pub n_iterations: usize,
}
/// A parsed OpenFOAM dictionary entry.
#[derive(Debug, Clone, PartialEq)]
pub enum FoamValue {
    /// A scalar number.
    Scalar(f64),
    /// A string or keyword.
    Word(String),
    /// A vector `(x y z)`.
    Vector([f64; 3]),
    /// A nested dictionary.
    Dict(FoamDict),
    /// A list of values.
    List(Vec<FoamValue>),
}
/// Writer for `constant/transportProperties`.
pub struct TransportProperties {
    /// Transport model: "Newtonian", "CrossPowerLaw", etc.
    pub transport_model: String,
    /// Kinematic viscosity \[m^2/s\].
    pub nu: f64,
}
impl TransportProperties {
    /// Create Newtonian transport with given kinematic viscosity.
    pub fn newtonian(nu: f64) -> Self {
        TransportProperties {
            transport_model: "Newtonian".to_string(),
            nu,
        }
    }
}
impl std::fmt::Display for TransportProperties {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = foam_header("dictionary", "transportProperties");
        s.push('\n');
        s.push_str(&format!("transportModel  {};\n\n", self.transport_model));
        s.push_str(&format!("nu              [0 2 -1 0 0 0 0] {};\n", self.nu));
        write!(f, "{s}")
    }
}
/// A scalar or vector field for OpenFOAM output (volScalarField / volVectorField).
pub struct FoamField {
    /// Number of internal cells.
    pub n_cells: usize,
    /// Field name used in the file header (e.g. "p", "U").
    pub field_name: String,
    /// OpenFOAM class: "volScalarField" or "volVectorField".
    pub field_class: String,
    /// Dimension set string, e.g. `"[0 2 -2 0 0 0 0]"` for pressure.
    pub dimensions: String,
    /// Internal field values.
    pub internal_values: FieldValues,
    /// Boundary condition list (one entry per patch).
    pub boundary_conditions: Vec<FoamBc>,
}
impl FoamField {}
impl std::fmt::Display for FoamField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = foam_header(&self.field_class, &self.field_name);
        s.push('\n');
        s.push_str(&format!("dimensions      {};\n\n", self.dimensions));
        match &self.internal_values {
            FieldValues::Uniform(v) => {
                s.push_str(&format!("internalField   uniform {};\n\n", v));
            }
            FieldValues::UniformVec(v) => {
                s.push_str(&format!(
                    "internalField   uniform ({} {} {});\n\n",
                    v[0], v[1], v[2]
                ));
            }
            FieldValues::NonUniform(vals) => {
                s.push_str("internalField   nonuniform List<scalar>\n");
                s.push_str(&format!("{}\n(\n", vals.len()));
                for v in vals {
                    s.push_str(&format!("{}\n", v));
                }
                s.push_str(");\n\n");
            }
            FieldValues::NonUniformVec(vals) => {
                s.push_str("internalField   nonuniform List<vector>\n");
                s.push_str(&format!("{}\n(\n", vals.len()));
                for v in vals {
                    s.push_str(&format!("({} {} {})\n", v[0], v[1], v[2]));
                }
                s.push_str(");\n\n");
            }
        }
        s.push_str("boundaryField\n{\n");
        for bc in &self.boundary_conditions {
            s.push_str(&format!("    {}\n    {{\n", bc.patch_name));
            s.push_str(&format!("        type        {};\n", bc.bc_type));
            if let Some(val) = &bc.value {
                s.push_str(&format!("        value       {};\n", val));
            }
            s.push_str("    }\n");
        }
        s.push_str("}\n");
        write!(f, "{s}")
    }
}
/// Internal/boundary field values for a `FoamField`.
pub enum FieldValues {
    /// Spatially uniform scalar value.
    Uniform(f64),
    /// Spatially uniform vector value.
    UniformVec([f64; 3]),
    /// Non-uniform scalar field (one value per cell).
    NonUniform(Vec<f64>),
    /// Non-uniform vector field (one 3-component vector per cell).
    NonUniformVec(Vec<[f64; 3]>),
}
/// Writer for the OpenFOAM `system/fvSolution` file.
pub struct FvSolution {
    /// Solver settings per field (field_name, solver, preconditioner, tolerance, relTol).
    pub solvers: Vec<FvSolverEntry>,
    /// SIMPLE/PISO/PIMPLE algorithm settings.
    pub algorithm: String,
    /// Number of correctors for PISO/PIMPLE.
    pub n_correctors: usize,
    /// Number of non-orthogonal correctors.
    pub n_non_orthogonal_correctors: usize,
    /// Pressure reference cell.
    pub p_ref_cell: usize,
    /// Pressure reference value.
    pub p_ref_value: f64,
}
impl FvSolution {
    /// Create default PISO solution settings.
    pub fn default_piso() -> Self {
        FvSolution {
            solvers: vec![
                FvSolverEntry {
                    field_name: "p".to_string(),
                    solver: "PCG".to_string(),
                    preconditioner: "DIC".to_string(),
                    tolerance: 1e-6,
                    rel_tol: 0.05,
                },
                FvSolverEntry {
                    field_name: "U".to_string(),
                    solver: "smoothSolver".to_string(),
                    preconditioner: "symGaussSeidel".to_string(),
                    tolerance: 1e-5,
                    rel_tol: 0.0,
                },
            ],
            algorithm: "PISO".to_string(),
            n_correctors: 2,
            n_non_orthogonal_correctors: 0,
            p_ref_cell: 0,
            p_ref_value: 0.0,
        }
    }
}
impl std::fmt::Display for FvSolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = foam_header("dictionary", "fvSolution");
        s.push('\n');
        s.push_str("solvers\n{\n");
        for entry in &self.solvers {
            s.push_str(&format!("    {}\n    {{\n", entry.field_name));
            s.push_str(&format!("        solver          {};\n", entry.solver));
            s.push_str(&format!(
                "        preconditioner  {};\n",
                entry.preconditioner
            ));
            s.push_str(&format!("        tolerance       {};\n", entry.tolerance));
            s.push_str(&format!("        relTol          {};\n", entry.rel_tol));
            s.push_str("    }\n");
        }
        s.push_str("}\n\n");
        s.push_str(&format!("{}\n{{\n", self.algorithm));
        s.push_str(&format!("    nCorrectors     {};\n", self.n_correctors));
        s.push_str(&format!(
            "    nNonOrthogonalCorrectors {};\n",
            self.n_non_orthogonal_correctors
        ));
        s.push_str(&format!("    pRefCell        {};\n", self.p_ref_cell));
        s.push_str(&format!("    pRefValue       {};\n", self.p_ref_value));
        s.push_str("}\n");
        write!(f, "{s}")
    }
}
/// A single solver entry in fvSolution.
#[derive(Debug, Clone)]
pub struct FvSolverEntry {
    /// Field name (e.g. "p", "U").
    pub field_name: String,
    /// Linear solver name (e.g. "PCG", "smoothSolver").
    pub solver: String,
    /// Preconditioner name (e.g. "DIC", "symGaussSeidel").
    pub preconditioner: String,
    /// Absolute convergence tolerance.
    pub tolerance: f64,
    /// Relative tolerance.
    pub rel_tol: f64,
}
/// Parsed contents of an OpenFOAM `FoamFile` header block.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamFileHeader {
    /// Format version (e.g., `2.0`).
    pub version: f64,
    /// Format string: `"ascii"` or `"binary"`.
    pub format: String,
    /// OpenFOAM class (e.g., `"volScalarField"`, `"polyMesh"`).
    pub class: String,
    /// Object name (e.g., `"p"`, `"points"`).
    pub object: String,
    /// Optional note/description field.
    pub note: Option<String>,
    /// Optional location path.
    pub location: Option<String>,
}
impl FoamFileHeader {
    /// Parse a `FoamFile { ... }` header from OpenFOAM file content.
    ///
    /// Returns `None` if no `FoamFile` block is found.
    pub fn parse(input: &str) -> Option<Self> {
        let cleaned = strip_foam_comments(input);
        let start = cleaned.find("FoamFile")?;
        let after = &cleaned[start + "FoamFile".len()..];
        let brace_start = after.find('{')?;
        let brace_content = &after[brace_start + 1..];
        let brace_end = brace_content.find('}')?;
        let block = &brace_content[..brace_end];
        let mut version = 2.0_f64;
        let mut format = "ascii".to_string();
        let mut class = String::new();
        let mut object = String::new();
        let mut note = None;
        let mut location = None;
        for line in block.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with("//") {
                continue;
            }
            let line_clean = trimmed.trim_end_matches(';').trim();
            let parts: Vec<&str> = line_clean.splitn(2, char::is_whitespace).collect();
            if parts.len() < 2 {
                continue;
            }
            let key = parts[0].trim();
            let val = parts[1].trim().trim_matches('"');
            match key {
                "version" => {
                    version = val.parse().unwrap_or(2.0);
                }
                "format" => {
                    format = val.to_string();
                }
                "class" => {
                    class = val.to_string();
                }
                "object" => {
                    object = val.to_string();
                }
                "note" => {
                    note = Some(val.to_string());
                }
                "location" => {
                    location = Some(val.to_string());
                }
                _ => {}
            }
        }
        if class.is_empty() && object.is_empty() {
            return None;
        }
        Some(FoamFileHeader {
            version,
            format,
            class,
            object,
            note,
            location,
        })
    }
}
impl std::fmt::Display for FoamFileHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        s.push_str("FoamFile\n{\n");
        s.push_str(&format!("    version     {};\n", self.version));
        s.push_str(&format!("    format      {};\n", self.format));
        s.push_str(&format!("    class       {};\n", self.class));
        if let Some(ref loc) = self.location {
            s.push_str(&format!("    location    \"{}\";\n", loc));
        }
        s.push_str(&format!("    object      {};\n", self.object));
        if let Some(ref note) = self.note {
            s.push_str(&format!("    note        \"{}\";\n", note));
        }
        s.push_str("}\n");
        write!(f, "{s}")
    }
}
