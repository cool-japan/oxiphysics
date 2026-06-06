//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

/// A set of anatomical landmarks.
#[derive(Clone, Debug)]
pub struct LandmarkSet {
    /// Collection of landmarks.
    pub landmarks: Vec<Landmark>,
    /// Coordinate space name.
    pub space: String,
}
impl LandmarkSet {
    /// Create an empty landmark set.
    pub fn new(space: &str) -> Self {
        Self {
            landmarks: Vec::new(),
            space: space.to_string(),
        }
    }
    /// Add a landmark.
    pub fn add(&mut self, lm: Landmark) {
        self.landmarks.push(lm);
    }
    /// Find landmark by name.
    pub fn find(&self, name: &str) -> Option<&Landmark> {
        self.landmarks.iter().find(|l| l.name == name)
    }
    /// Compute centroid of all landmarks.
    pub fn centroid(&self) -> Option<[f64; 3]> {
        if self.landmarks.is_empty() {
            return None;
        }
        let n = self.landmarks.len() as f64;
        let mut c = [0.0f64; 3];
        for lm in &self.landmarks {
            for (cd, &p) in c.iter_mut().zip(lm.position.iter()) {
                *cd += p;
            }
        }
        for cd in c.iter_mut() {
            *cd /= n;
        }
        Some(c)
    }
}
/// DICOM data element: tag + typed value.
#[derive(Clone, Debug)]
pub struct DicomElement {
    /// The tag identifying this element.
    pub tag: DicomTag,
    /// The value representation and data.
    pub vr: DicomVr,
}
impl DicomElement {
    /// Create a new DICOM element.
    pub fn new(tag: DicomTag, vr: DicomVr) -> Self {
        Self { tag, vr }
    }
}
/// VTK structured image data (vtkImageData equivalent).
#[derive(Clone, Debug)]
pub struct VtkImageData {
    /// Grid dimensions \[nx, ny, nz\].
    pub dimensions: [usize; 3],
    /// Grid origin \[ox, oy, oz\].
    pub origin: [f64; 3],
    /// Voxel spacing \[dx, dy, dz\].
    pub spacing: [f64; 3],
    /// Named scalar arrays (one value per voxel).
    pub scalar_arrays: HashMap<String, Vec<f64>>,
    /// Named vector arrays (3 values per voxel).
    pub vector_arrays: HashMap<String, Vec<[f64; 3]>>,
}
impl VtkImageData {
    /// Create new VTK image data.
    pub fn new(dimensions: [usize; 3], origin: [f64; 3], spacing: [f64; 3]) -> Self {
        Self {
            dimensions,
            origin,
            spacing,
            scalar_arrays: HashMap::new(),
            vector_arrays: HashMap::new(),
        }
    }
    /// Total number of voxels.
    pub fn n_voxels(&self) -> usize {
        self.dimensions[0] * self.dimensions[1] * self.dimensions[2]
    }
    /// Add a scalar data array.
    pub fn add_scalar(&mut self, name: &str, data: Vec<f64>) {
        self.scalar_arrays.insert(name.to_string(), data);
    }
    /// Add a vector data array.
    pub fn add_vector(&mut self, name: &str, data: Vec<[f64; 3]>) {
        self.vector_arrays.insert(name.to_string(), data);
    }
    /// Convert voxel index to world coordinates.
    pub fn index_to_world(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + i as f64 * self.spacing[0],
            self.origin[1] + j as f64 * self.spacing[1],
            self.origin[2] + k as f64 * self.spacing[2],
        ]
    }
    /// Serialize to VTK legacy format string.
    pub fn to_vtk_string(&self, scalar_name: &str) -> String {
        let [nx, ny, nz] = self.dimensions;
        let [ox, oy, oz] = self.origin;
        let [dx, dy, dz] = self.spacing;
        let n = self.n_voxels();
        let mut s = format!(
            "# vtk DataFile Version 3.0\nVTK ImageData\nASCII\nDATASET STRUCTURED_POINTS\n\
             DIMENSIONS {} {} {}\nORIGIN {} {} {}\nSPACING {} {} {}\n\
             POINT_DATA {}\n",
            nx, ny, nz, ox, oy, oz, dx, dy, dz, n
        );
        if let Some(arr) = self.scalar_arrays.get(scalar_name) {
            s.push_str(&format!(
                "SCALARS {} float 1\nLOOKUP_TABLE default\n",
                scalar_name
            ));
            for &v in arr.iter().take(10) {
                s.push_str(&format!("{:.6} ", v));
            }
        }
        s
    }
}
/// NIfTI qform code.
#[derive(Clone, Debug, PartialEq)]
pub enum NiftiQformCode {
    /// Unknown / not set.
    Unknown,
    /// Scanner anatomical coordinates.
    ScannerAnat,
    /// Aligned to reference.
    AlignedAnat,
    /// Talairach space.
    Talairach,
    /// MNI 152 standard space.
    Mni152,
}
impl NiftiQformCode {
    /// Return the integer code as stored in the NIfTI header byte.
    pub fn as_code(&self) -> u16 {
        match self {
            NiftiQformCode::Unknown => 0,
            NiftiQformCode::ScannerAnat => 1,
            NiftiQformCode::AlignedAnat => 2,
            NiftiQformCode::Talairach => 3,
            NiftiQformCode::Mni152 => 4,
        }
    }
    /// Create from integer code.
    pub fn from_code(code: u16) -> Self {
        match code {
            1 => NiftiQformCode::ScannerAnat,
            2 => NiftiQformCode::AlignedAnat,
            3 => NiftiQformCode::Talairach,
            4 => NiftiQformCode::Mni152,
            _ => NiftiQformCode::Unknown,
        }
    }
}
/// NRRD (Nearly Raw Raster Data) data encoding.
#[derive(Clone, Debug, PartialEq)]
pub enum NrrdEncoding {
    /// Raw binary (no compression).
    Raw,
    /// ASCII text.
    Text,
    /// Gzip compressed (mock: same as raw in unit tests).
    Gzip,
}
/// DICOM header tag (group, element).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DicomTag {
    /// Tag group number.
    pub group: u16,
    /// Tag element number.
    pub element: u16,
}
impl DicomTag {
    /// Create a new DICOM tag.
    pub fn new(group: u16, element: u16) -> Self {
        Self { group, element }
    }
    /// Patient name tag (0010,0010).
    pub fn patient_name() -> Self {
        Self::new(0x0010, 0x0010)
    }
    /// Patient ID tag (0010,0020).
    pub fn patient_id() -> Self {
        Self::new(0x0010, 0x0020)
    }
    /// Modality tag (0008,0060).
    pub fn modality() -> Self {
        Self::new(0x0008, 0x0060)
    }
    /// Rows tag (0028,0010).
    pub fn rows() -> Self {
        Self::new(0x0028, 0x0010)
    }
    /// Columns tag (0028,0011).
    pub fn columns() -> Self {
        Self::new(0x0028, 0x0011)
    }
    /// Pixel spacing tag (0028,0030).
    pub fn pixel_spacing() -> Self {
        Self::new(0x0028, 0x0030)
    }
    /// Slice thickness tag (0050,0018).
    pub fn slice_thickness() -> Self {
        Self::new(0x0050, 0x0018)
    }
    /// Window center tag (0028,1050).
    pub fn window_center() -> Self {
        Self::new(0x0028, 0x1050)
    }
    /// Window width tag (0028,1051).
    pub fn window_width() -> Self {
        Self::new(0x0028, 0x1051)
    }
}
/// NIfTI-1 header (.nii / .nii.gz reader mock).
#[derive(Clone, Debug)]
pub struct NiftiReader {
    /// File path.
    pub file_path: String,
    /// Image dimensions \[nx, ny, nz, nt, ...\].
    pub dim: [usize; 8],
    /// Voxel sizes in mm.
    pub pixdim: [f64; 8],
    /// Data type.
    pub datatype: NiftiDtype,
    /// Affine transform (4x4, row-major).
    pub affine: [f64; 16],
    /// Description string.
    pub descrip: String,
    /// Voxel data (as f64 for convenience).
    pub data: Vec<f64>,
}
impl NiftiReader {
    /// Create a new NIfTI reader.
    pub fn new(file_path: &str) -> Self {
        let mut affine = [0.0f64; 16];
        affine[0] = 1.0;
        affine[5] = 1.0;
        affine[10] = 1.0;
        affine[15] = 1.0;
        Self {
            file_path: file_path.to_string(),
            dim: [3, 64, 64, 32, 1, 1, 1, 1],
            pixdim: [1.0; 8],
            datatype: NiftiDtype::Float32,
            affine,
            descrip: String::new(),
            data: Vec::new(),
        }
    }
    /// Number of voxels = dim\[1\] * dim\[2\] * dim\[3\].
    pub fn n_voxels(&self) -> usize {
        self.dim[1] * self.dim[2] * self.dim[3]
    }
    /// Initialize data with zeros.
    pub fn init_data(&mut self) {
        self.data = vec![0.0; self.n_voxels()];
    }
    /// Apply affine to voxel index (i, j, k) to get world coordinates (mm).
    pub fn voxel_to_world(&self, i: f64, j: f64, k: f64) -> [f64; 3] {
        let a = &self.affine;
        [
            a[0] * i + a[1] * j + a[2] * k + a[3],
            a[4] * i + a[5] * j + a[6] * k + a[7],
            a[8] * i + a[9] * j + a[10] * k + a[11],
        ]
    }
    /// Get voxel value at (i, j, k).
    pub fn get_voxel(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        let nx = self.dim[1];
        let ny = self.dim[2];
        let nz = self.dim[3];
        if i < nx && j < ny && k < nz {
            let idx = k * nx * ny + j * nx + i;
            self.data.get(idx).copied()
        } else {
            None
        }
    }
}
/// An anatomical landmark.
#[derive(Clone, Debug)]
pub struct Landmark {
    /// Landmark name (e.g., "AC", "PC").
    pub name: String,
    /// World coordinates in mm.
    pub position: [f64; 3],
    /// Uncertainty radius (mm).
    pub uncertainty: f64,
    /// Whether this landmark was manually placed.
    pub manual: bool,
}
impl Landmark {
    /// Create a new landmark.
    pub fn new(name: &str, position: [f64; 3], uncertainty: f64) -> Self {
        Self {
            name: name.to_string(),
            position,
            uncertainty,
            manual: true,
        }
    }
}
/// NIfTI sform/qform affine parameters.
#[derive(Clone, Debug)]
pub struct NiftiTransform {
    /// qform code.
    pub qform_code: NiftiQformCode,
    /// sform code.
    pub sform_code: u16,
    /// Quaternion parameters (qb, qc, qd) for qform.
    pub quatern: [f64; 3],
    /// qoffset (qx, qy, qz).
    pub qoffset: [f64; 3],
    /// Pixel dimension sign (qfac: +1 or -1).
    pub qfac: f64,
    /// sform row vectors for the 3×4 affine (srow_x, srow_y, srow_z).
    pub sform_matrix: [[f64; 4]; 3],
}
impl NiftiTransform {
    /// Create an identity-like transform.
    pub fn identity() -> Self {
        Self {
            qform_code: NiftiQformCode::Unknown,
            sform_code: 0,
            quatern: [0.0, 0.0, 0.0],
            qoffset: [0.0, 0.0, 0.0],
            qfac: 1.0,
            sform_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
        }
    }
    /// Convert voxel index (i,j,k) to world coordinates using the sform matrix.
    pub fn sform_to_world(&self, i: f64, j: f64, k: f64) -> [f64; 3] {
        let m = &self.sform_matrix;
        [
            m[0][0] * i + m[0][1] * j + m[0][2] * k + m[0][3],
            m[1][0] * i + m[1][1] * j + m[1][2] * k + m[1][3],
            m[2][0] * i + m[2][1] * j + m[2][2] * k + m[2][3],
        ]
    }
}
/// NIfTI-1 data types.
#[derive(Clone, Debug, PartialEq)]
pub enum NiftiDtype {
    /// Unsigned 8-bit integer.
    Uint8,
    /// Signed 16-bit integer.
    Int16,
    /// Signed 32-bit integer.
    Int32,
    /// 32-bit float.
    Float32,
    /// 64-bit float.
    Float64,
}
/// A triangle in 3D (for STL output).
#[derive(Clone, Debug)]
pub struct StlTriangle {
    /// Triangle normal (unit vector).
    pub normal: [f32; 3],
    /// Three vertex positions.
    pub vertices: [[f32; 3]; 3],
}
impl StlTriangle {
    /// Compute the face normal from vertex positions.
    pub fn compute_normal(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> [f32; 3] {
        let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-30);
        [n[0] / len, n[1] / len, n[2] / len]
    }
}
/// Triangulated surface mesh from medical images.
#[derive(Clone, Debug)]
pub struct MedicalMesh {
    /// Vertex positions (world coordinates, mm).
    pub vertices: Vec<[f64; 3]>,
    /// Triangle connectivity (indices into vertices).
    pub triangles: Vec<[usize; 3]>,
    /// Vertex normals.
    pub normals: Vec<[f64; 3]>,
    /// Mesh name / label.
    pub name: String,
}
impl MedicalMesh {
    /// Create an empty mesh.
    pub fn new(name: &str) -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
            normals: Vec::new(),
            name: name.to_string(),
        }
    }
    /// Add a triangle.
    pub fn add_triangle(&mut self, v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) {
        let i0 = self.vertices.len();
        self.vertices.push(v0);
        self.vertices.push(v1);
        self.vertices.push(v2);
        self.triangles.push([i0, i0 + 1, i0 + 2]);
        let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-30);
        let nn = [n[0] / len, n[1] / len, n[2] / len];
        self.normals.push(nn);
        self.normals.push(nn);
        self.normals.push(nn);
    }
    /// Compute surface area.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0;
        for &[i0, i1, i2] in &self.triangles {
            let v0 = self.vertices[i0];
            let v1 = self.vertices[i1];
            let v2 = self.vertices[i2];
            let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let cx = e1[1] * e2[2] - e1[2] * e2[1];
            let cy = e1[2] * e2[0] - e1[0] * e2[2];
            let cz = e1[0] * e2[1] - e1[1] * e2[0];
            area += 0.5 * (cx * cx + cy * cy + cz * cz).sqrt();
        }
        area
    }
    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }
}
/// Reconstruct a 3D voxel volume from a stack of 2D slices.
///
/// Each slice is a `Vec`u16` of `rows × cols` pixels.  The function stacks
/// them into a single flat buffer in slice-major order (slice fastest varies last).
#[derive(Clone, Debug)]
pub struct VolumeReconstructor {
    /// Number of rows per slice.
    pub rows: usize,
    /// Number of columns per slice.
    pub cols: usize,
    /// Pixel spacing [row_mm, col_mm].
    pub pixel_spacing: [f64; 2],
    /// Slice thickness in mm.
    pub slice_thickness: f64,
    /// Assembled voxel data (slice-major).
    pub data: Vec<u16>,
    /// Number of slices assembled so far.
    pub n_slices: usize,
}
impl VolumeReconstructor {
    /// Create a new reconstructor for slices of size `rows × cols`.
    pub fn new(rows: usize, cols: usize, pixel_spacing: [f64; 2], slice_thickness: f64) -> Self {
        Self {
            rows,
            cols,
            pixel_spacing,
            slice_thickness,
            data: Vec::new(),
            n_slices: 0,
        }
    }
    /// Append one slice of pixel data.
    pub fn add_slice(&mut self, slice: &[u16]) {
        assert_eq!(slice.len(), self.rows * self.cols, "slice size mismatch");
        self.data.extend_from_slice(slice);
        self.n_slices += 1;
    }
    /// Access voxel at (row, col, slice).
    pub fn voxel(&self, row: usize, col: usize, slice: usize) -> Option<u16> {
        if row < self.rows && col < self.cols && slice < self.n_slices {
            Some(self.data[slice * self.rows * self.cols + row * self.cols + col])
        } else {
            None
        }
    }
    /// Physical size in mm [rows, cols, slices].
    pub fn physical_size(&self) -> [f64; 3] {
        [
            self.rows as f64 * self.pixel_spacing[0],
            self.cols as f64 * self.pixel_spacing[1],
            self.n_slices as f64 * self.slice_thickness,
        ]
    }
}
/// NRRD header and data container.
#[derive(Clone, Debug)]
pub struct NrrdReader {
    /// Key-value pairs from the NRRD header.
    pub fields: HashMap<String, String>,
    /// Data encoding.
    pub encoding: NrrdEncoding,
    /// Flat f64 data buffer.
    pub data: Vec<f64>,
    /// Dimensions as parsed from the "sizes" field.
    pub sizes: Vec<usize>,
}
impl NrrdReader {
    /// Create an empty NRRD reader.
    pub fn new() -> Self {
        Self {
            fields: HashMap::new(),
            encoding: NrrdEncoding::Raw,
            data: Vec::new(),
            sizes: Vec::new(),
        }
    }
    /// Parse an NRRD header string (lines of "key: value").
    ///
    /// Populates `fields`, `encoding`, and `sizes`.
    pub fn parse_header(&mut self, header: &str) {
        for line in header.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some((k, v)) = line.split_once(':') {
                let key = k.trim().to_lowercase();
                let val = v.trim().to_string();
                if key == "encoding" {
                    self.encoding = match val.as_str() {
                        "raw" => NrrdEncoding::Raw,
                        "text" | "ascii" => NrrdEncoding::Text,
                        "gzip" | "gz" => NrrdEncoding::Gzip,
                        _ => NrrdEncoding::Raw,
                    };
                }
                if key == "sizes" {
                    self.sizes = val
                        .split_whitespace()
                        .filter_map(|s| s.parse().ok())
                        .collect();
                }
                self.fields.insert(key, val);
            }
        }
    }
    /// Load raw-binary f32 data from a byte slice.
    pub fn load_raw_f32(&mut self, bytes: &[u8]) {
        self.data = bytes
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64)
            .collect();
    }
    /// Load text-encoded data (space/newline-separated ASCII floats).
    pub fn load_text(&mut self, text: &str) {
        self.data = text
            .split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();
    }
    /// Total number of voxels as the product of `sizes`.
    pub fn n_voxels(&self) -> usize {
        self.sizes.iter().product::<usize>()
    }
}
/// MINC volume (mock in-memory representation).
#[derive(Clone, Debug)]
pub struct MincVolume {
    /// Spatial dimensions (up to 4: x, y, z, t).
    pub dimensions: Vec<MincDimension>,
    /// Flat data buffer.
    pub data: Vec<f64>,
    /// Real-value valid range [min, max].
    pub valid_range: [f64; 2],
    /// Attribute map (e.g. history, patient info).
    pub attributes: HashMap<String, String>,
}
impl MincVolume {
    /// Create a new MINC volume.
    pub fn new(dimensions: Vec<MincDimension>) -> Self {
        let n: usize = dimensions.iter().map(|d| d.length).product();
        Self {
            dimensions,
            data: vec![0.0; n],
            valid_range: [0.0, 1.0],
            attributes: HashMap::new(),
        }
    }
    /// Total voxel count.
    pub fn n_voxels(&self) -> usize {
        self.dimensions.iter().map(|d| d.length).product()
    }
    /// Set attribute.
    pub fn set_attr(&mut self, key: &str, value: &str) {
        self.attributes.insert(key.to_string(), value.to_string());
    }
    /// Voxel index for (i, j, k) in a 3-D volume (last-index fastest = C order).
    pub fn index3(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        if self.dimensions.len() < 3 {
            return None;
        }
        let nx = self.dimensions[0].length;
        let ny = self.dimensions[1].length;
        let nz = self.dimensions[2].length;
        if i < nx && j < ny && k < nz {
            Some(i * ny * nz + j * nz + k)
        } else {
            None
        }
    }
}
/// MetaImage (.mha / .mhd) format.
#[derive(Clone, Debug)]
pub struct MhaMhdFormat {
    /// Header key-value pairs.
    pub header: HashMap<String, String>,
    /// Raw voxel data.
    pub data: Vec<f64>,
    /// Whether the data is embedded (.mha) or external (.mhd).
    pub is_mha: bool,
}
impl MhaMhdFormat {
    /// Create a new MetaImage container.
    pub fn new(is_mha: bool) -> Self {
        let mut header = HashMap::new();
        header.insert("ObjectType".to_string(), "Image".to_string());
        header.insert("NDims".to_string(), "3".to_string());
        header.insert("ElementType".to_string(), "MET_FLOAT".to_string());
        Self {
            header,
            data: Vec::new(),
            is_mha,
        }
    }
    /// Set header field.
    pub fn set(&mut self, key: &str, value: &str) {
        self.header.insert(key.to_string(), value.to_string());
    }
    /// Get header field.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.header.get(key).map(|s| s.as_str())
    }
    /// Parse dimensions from "DimSize" header field.
    pub fn dimensions(&self) -> Option<Vec<usize>> {
        self.get("DimSize").map(|s| {
            s.split_whitespace()
                .filter_map(|t| t.parse().ok())
                .collect()
        })
    }
    /// Serialize header to string.
    pub fn serialize_header(&self) -> String {
        let mut lines: Vec<String> = self
            .header
            .iter()
            .map(|(k, v)| format!("{} = {}", k, v))
            .collect();
        lines.sort();
        if self.is_mha {
            lines.push("ElementDataFile = LOCAL".to_string());
        }
        lines.join("\n")
    }
}
/// DICOM header reader (mock implementation).
#[derive(Clone, Debug)]
pub struct DicomReader {
    /// File path (virtual).
    pub file_path: String,
    /// Header key-value pairs (string representation).
    pub header: HashMap<DicomTag, String>,
    /// Pixel data as u16 values (row-major).
    pub pixel_data: Vec<u16>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub columns: usize,
}
impl DicomReader {
    /// Create a new DICOM reader.
    pub fn new(file_path: &str) -> Self {
        Self {
            file_path: file_path.to_string(),
            header: HashMap::new(),
            pixel_data: Vec::new(),
            rows: 0,
            columns: 0,
        }
    }
    /// Set a header tag value.
    pub fn set_tag(&mut self, tag: DicomTag, value: &str) {
        self.header.insert(tag, value.to_string());
    }
    /// Get a header tag value.
    pub fn get_tag(&self, tag: &DicomTag) -> Option<&str> {
        self.header.get(tag).map(|s| s.as_str())
    }
    /// Load synthetic pixel data for testing (checker pattern).
    pub fn load_synthetic(&mut self, rows: usize, cols: usize) {
        self.rows = rows;
        self.columns = cols;
        self.pixel_data = (0..(rows * cols))
            .map(|i| {
                let r = i / cols;
                let c = i % cols;
                if (r + c).is_multiple_of(2) {
                    2048u16
                } else {
                    512u16
                }
            })
            .collect();
        self.set_tag(DicomTag::rows(), &rows.to_string());
        self.set_tag(DicomTag::columns(), &cols.to_string());
    }
    /// Extract pixel at (row, col).
    pub fn pixel(&self, row: usize, col: usize) -> Option<u16> {
        if row < self.rows && col < self.columns {
            Some(self.pixel_data[row * self.columns + col])
        } else {
            None
        }
    }
    /// Get patient name from header.
    pub fn patient_name(&self) -> Option<&str> {
        self.get_tag(&DicomTag::patient_name())
    }
    /// Get modality from header.
    pub fn modality(&self) -> Option<&str> {
        self.get_tag(&DicomTag::modality())
    }
}
/// NIfTI-1 writer.
#[derive(Clone, Debug)]
pub struct NiftiWriter {
    /// Output path.
    pub output_path: String,
    /// Data scaling: slope.
    pub scl_slope: f64,
    /// Data scaling: intercept.
    pub scl_inter: f64,
    /// Use gzip compression.
    pub gzip: bool,
}
impl NiftiWriter {
    /// Create a new NIfTI writer.
    pub fn new(output_path: &str) -> Self {
        Self {
            output_path: output_path.to_string(),
            scl_slope: 1.0,
            scl_inter: 0.0,
            gzip: false,
        }
    }
    /// Build NIfTI-1 header bytes (simplified, 348-byte mock).
    pub fn build_header(&self, reader: &NiftiReader) -> Vec<u8> {
        let mut header = vec![0u8; 348];
        let sizeof_hdr = 348u32.to_le_bytes();
        header[0] = sizeof_hdr[0];
        header[1] = sizeof_hdr[1];
        let ndim = reader.dim[0] as u16;
        let bytes = ndim.to_le_bytes();
        header[40] = bytes[0];
        header[41] = bytes[1];
        header[344] = b'n';
        header[345] = b'+';
        header[346] = b'1';
        header[347] = 0;
        header
    }
    /// Write NIfTI data to a buffer (mock: returns serialized bytes).
    pub fn write_to_buffer(&self, reader: &NiftiReader) -> Vec<u8> {
        let header = self.build_header(reader);
        let mut buf = header;
        for &v in &reader.data {
            let scaled = (v * self.scl_slope + self.scl_inter) as f32;
            buf.extend_from_slice(&scaled.to_le_bytes());
        }
        buf
    }
}
/// A binary or multi-label 3D segmentation mask.
#[derive(Clone, Debug)]
pub struct SegmentationMask {
    /// Label array (0 = background).
    pub labels: Vec<u8>,
    /// Dimensions [nx, ny, nz].
    pub dimensions: [usize; 3],
    /// Voxel spacing [dx, dy, dz] in mm.
    pub spacing: [f64; 3],
    /// Label names.
    pub label_names: HashMap<u8, String>,
}
impl SegmentationMask {
    /// Create a new empty segmentation mask.
    pub fn new(dimensions: [usize; 3], spacing: [f64; 3]) -> Self {
        let n = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            labels: vec![0; n],
            dimensions,
            spacing,
            label_names: HashMap::new(),
        }
    }
    /// Set a voxel label.
    pub fn set_label(&mut self, i: usize, j: usize, k: usize, label: u8) {
        let [nx, ny, _nz] = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        if idx < self.labels.len() {
            self.labels[idx] = label;
        }
    }
    /// Get a voxel label.
    pub fn get_label(&self, i: usize, j: usize, k: usize) -> Option<u8> {
        let [nx, ny, _nz] = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        self.labels.get(idx).copied()
    }
    /// Count voxels with a given label.
    pub fn count_label(&self, label: u8) -> usize {
        self.labels.iter().filter(|&&l| l == label).count()
    }
    /// Volume of labeled region in mm^3.
    pub fn volume_mm3(&self, label: u8) -> f64 {
        let n = self.count_label(label) as f64;
        let voxel_vol = self.spacing[0] * self.spacing[1] * self.spacing[2];
        n * voxel_vol
    }
    /// Bounding box of labeled region [min_i, max_i, min_j, max_j, min_k, max_k].
    pub fn bounding_box(&self, label: u8) -> Option<[usize; 6]> {
        let [nx, ny, nz] = self.dimensions;
        let mut min_i = nx;
        let mut max_i = 0;
        let mut min_j = ny;
        let mut max_j = 0;
        let mut min_k = nz;
        let mut max_k = 0;
        let mut found = false;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if self.get_label(i, j, k) == Some(label) {
                        min_i = min_i.min(i);
                        max_i = max_i.max(i);
                        min_j = min_j.min(j);
                        max_j = max_j.max(j);
                        min_k = min_k.min(k);
                        max_k = max_k.max(k);
                        found = true;
                    }
                }
            }
        }
        if found {
            Some([min_i, max_i, min_j, max_j, min_k, max_k])
        } else {
            None
        }
    }
    /// Register a label name.
    pub fn name_label(&mut self, label: u8, name: &str) {
        self.label_names.insert(label, name.to_string());
    }
}
/// DICOM data set: an ordered collection of elements.
#[derive(Clone, Debug, Default)]
pub struct DicomDataSet {
    /// Elements stored in tag order.
    pub elements: Vec<DicomElement>,
}
impl DicomDataSet {
    /// Create an empty data set.
    pub fn new() -> Self {
        Self {
            elements: Vec::new(),
        }
    }
    /// Insert or replace an element.
    pub fn set(&mut self, element: DicomElement) {
        if let Some(pos) = self.elements.iter().position(|e| e.tag == element.tag) {
            self.elements[pos] = element;
        } else {
            self.elements.push(element);
            self.elements.sort_by(|a, b| {
                a.tag
                    .group
                    .cmp(&b.tag.group)
                    .then(a.tag.element.cmp(&b.tag.element))
            });
        }
    }
    /// Look up an element by tag.
    pub fn get(&self, tag: &DicomTag) -> Option<&DicomElement> {
        self.elements.iter().find(|e| &e.tag == tag)
    }
    /// Convenience: get US value for a tag.
    pub fn get_us(&self, tag: &DicomTag) -> Option<u16> {
        self.get(tag)?.vr.as_us()
    }
    /// Convenience: get DS value for a tag.
    pub fn get_ds(&self, tag: &DicomTag) -> Option<f64> {
        self.get(tag)?.vr.as_ds()
    }
    /// Convenience: get string for LO/UI tag.
    pub fn get_str(&self, tag: &DicomTag) -> Option<&str> {
        self.get(tag)?.vr.as_str()
    }
    /// Parse a minimal DICOM-like byte stream (explicit VR, little-endian mock).
    ///
    /// Format: repeated records of `\[group:u16\]\[element:u16\]\[vr:2 bytes\]\[length:u16\]\[data\]`.
    /// This is a simplified subset suitable for unit-test round-trips.
    pub fn parse_bytes(bytes: &[u8]) -> Self {
        let mut ds = DicomDataSet::new();
        let mut pos = 0usize;
        while pos + 8 <= bytes.len() {
            let group = u16::from_le_bytes([bytes[pos], bytes[pos + 1]]);
            let element = u16::from_le_bytes([bytes[pos + 2], bytes[pos + 3]]);
            let vr_code = [bytes[pos + 4], bytes[pos + 5]];
            let length = u16::from_le_bytes([bytes[pos + 6], bytes[pos + 7]]) as usize;
            pos += 8;
            if pos + length > bytes.len() {
                break;
            }
            let data = &bytes[pos..pos + length];
            pos += length;
            let tag = DicomTag::new(group, element);
            let vr = match &vr_code {
                b"US" if length == 2 => DicomVr::Us(u16::from_le_bytes([data[0], data[1]])),
                b"DS" => {
                    let s = std::str::from_utf8(data).unwrap_or("0").trim();
                    DicomVr::Ds(s.parse::<f64>().unwrap_or(0.0))
                }
                b"LO" => DicomVr::Lo(String::from_utf8_lossy(data).trim().to_string()),
                b"UI" => DicomVr::Ui(String::from_utf8_lossy(data).trim().to_string()),
                b"OW" => DicomVr::OW(data.to_vec()),
                _ => continue,
            };
            ds.set(DicomElement::new(tag, vr));
        }
        ds
    }
    /// Serialize the data set back to bytes (same simplified format as `parse_bytes`).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for elem in &self.elements {
            let vr_code: &[u8; 2];
            let data: Vec<u8>;
            match &elem.vr {
                DicomVr::Us(v) => {
                    vr_code = b"US";
                    data = v.to_le_bytes().to_vec();
                }
                DicomVr::Ds(v) => {
                    vr_code = b"DS";
                    data = format!("{:.6}", v).into_bytes();
                }
                DicomVr::Lo(s) => {
                    vr_code = b"LO";
                    data = s.as_bytes().to_vec();
                }
                DicomVr::Ui(s) => {
                    vr_code = b"UI";
                    data = s.as_bytes().to_vec();
                }
                DicomVr::OW(bytes) => {
                    vr_code = b"OW";
                    data = bytes.clone();
                }
                DicomVr::Sq(_) => continue,
            }
            let length = data.len() as u16;
            out.extend_from_slice(&elem.tag.group.to_le_bytes());
            out.extend_from_slice(&elem.tag.element.to_le_bytes());
            out.extend_from_slice(vr_code);
            out.extend_from_slice(&length.to_le_bytes());
            out.extend_from_slice(&data);
        }
        out
    }
}
/// A single tractography fiber (a polyline of 3D points).
#[derive(Clone, Debug)]
pub struct FiberTract {
    /// Ordered points along the fiber [x, y, z] in mm.
    pub points: Vec<[f64; 3]>,
    /// Optional per-point scalar (e.g. FA, MD).
    pub scalars: Vec<f64>,
}
impl FiberTract {
    /// Create a new fiber tract.
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        let n = points.len();
        Self {
            points,
            scalars: vec![0.0; n],
        }
    }
    /// Arc length of the fiber in mm.
    pub fn arc_length(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        self.points
            .windows(2)
            .map(|w| {
                let d = [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .sum()
    }
    /// Number of points.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
}
/// Radiation therapy dose distribution.
#[derive(Clone, Debug)]
pub struct DoseVolume {
    /// Dose values in Gray (Gy).
    pub dose: Vec<f64>,
    /// Dimensions [nx, ny, nz].
    pub dimensions: [usize; 3],
    /// Voxel spacing [dx, dy, dz] in mm.
    pub spacing: [f64; 3],
    /// Prescription dose (Gy).
    pub prescription_dose: f64,
}
impl DoseVolume {
    /// Create a dose volume.
    pub fn new(dimensions: [usize; 3], spacing: [f64; 3], prescription_dose: f64) -> Self {
        let n = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            dose: vec![0.0; n],
            dimensions,
            spacing,
            prescription_dose,
        }
    }
    /// Get dose at voxel (i, j, k).
    pub fn get_dose(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        let [nx, ny, _nz] = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        self.dose.get(idx).copied()
    }
    /// Set dose at voxel.
    pub fn set_dose(&mut self, i: usize, j: usize, k: usize, d: f64) {
        let [nx, ny, _nz] = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        if idx < self.dose.len() {
            self.dose[idx] = d;
        }
    }
    /// Maximum dose in the volume.
    pub fn max_dose(&self) -> f64 {
        self.dose.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
    /// Mean dose over non-zero voxels.
    pub fn mean_dose(&self) -> f64 {
        let nonzero: Vec<f64> = self.dose.iter().cloned().filter(|&d| d > 0.0).collect();
        if nonzero.is_empty() {
            return 0.0;
        }
        nonzero.iter().sum::<f64>() / nonzero.len() as f64
    }
    /// Compute dose-volume histogram (DVH).
    ///
    /// Returns (dose_bins, volume_fraction) vectors.
    pub fn dvh(&self, n_bins: usize) -> (Vec<f64>, Vec<f64>) {
        let d_max = self.max_dose();
        if d_max <= 0.0 {
            return (vec![0.0; n_bins], vec![0.0; n_bins]);
        }
        let bin_width = d_max / n_bins as f64;
        let mut counts = vec![0usize; n_bins];
        let n_total = self.dose.len();
        for &d in &self.dose {
            let bin = ((d / d_max) * (n_bins - 1) as f64) as usize;
            counts[bin.min(n_bins - 1)] += 1;
        }
        let mut cumulative = vec![0usize; n_bins];
        let mut running = 0usize;
        for i in (0..n_bins).rev() {
            running += counts[i];
            cumulative[i] = running;
        }
        let dose_bins: Vec<f64> = (0..n_bins).map(|i| i as f64 * bin_width).collect();
        let vol_frac: Vec<f64> = cumulative
            .iter()
            .map(|&c| c as f64 / n_total as f64)
            .collect();
        (dose_bins, vol_frac)
    }
    /// V_n: volume fraction receiving at least `dose_threshold` Gy.
    pub fn v_dose(&self, dose_threshold: f64) -> f64 {
        let above = self.dose.iter().filter(|&&d| d >= dose_threshold).count();
        above as f64 / self.dose.len() as f64
    }
    /// D_n: dose received by at least `volume_fraction` of the volume.
    pub fn d_volume(&self, volume_fraction: f64) -> f64 {
        let mut sorted = self.dose.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idx = ((1.0 - volume_fraction) * sorted.len() as f64) as usize;
        sorted[idx.min(sorted.len() - 1)]
    }
}
/// VTK PolyData writer for tractography fibers.
pub struct TractVtkWriter {
    /// Label for the output data set.
    pub label: String,
}
impl TractVtkWriter {
    /// Create a new VTK tract writer.
    pub fn new(label: &str) -> Self {
        Self {
            label: label.to_string(),
        }
    }
    /// Serialise fibers to VTK legacy ASCII PolyData format.
    pub fn to_vtk_string(&self, fibers: &[FiberTract]) -> String {
        let total_pts: usize = fibers.iter().map(|f| f.n_points()).sum();
        let total_lines: usize = fibers.len();
        let cell_size: usize = fibers.iter().map(|f| f.n_points() + 1).sum();
        let mut s = format!(
            "# vtk DataFile Version 3.0\n{}\nASCII\nDATASET POLYDATA\nPOINTS {} float\n",
            self.label, total_pts
        );
        for fiber in fibers {
            for p in &fiber.points {
                s.push_str(&format!("{:.4} {:.4} {:.4}\n", p[0], p[1], p[2]));
            }
        }
        s.push_str(&format!(
            "LINES {} {}\n",
            total_lines,
            total_lines + cell_size
        ));
        let mut offset = 0usize;
        for fiber in fibers {
            s.push_str(&format!("{}", fiber.n_points()));
            for k in 0..fiber.n_points() {
                s.push_str(&format!(" {}", offset + k));
            }
            s.push('\n');
            offset += fiber.n_points();
        }
        if total_pts > 0 {
            s.push_str(&format!(
                "POINT_DATA {}\nSCALARS scalars float 1\nLOOKUP_TABLE default\n",
                total_pts
            ));
            for fiber in fibers {
                for &sc in &fiber.scalars {
                    s.push_str(&format!("{:.6} ", sc));
                }
            }
        }
        s
    }
}
/// DICOM image data with physical spacing and window/level settings.
#[derive(Clone, Debug)]
pub struct DicomImageData {
    /// Pixel array (row-major, in-plane first, slice last).
    pub pixels: Vec<u16>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Number of slices.
    pub n_slices: usize,
    /// Pixel spacing in mm [row_spacing, col_spacing].
    pub pixel_spacing: [f64; 2],
    /// Slice thickness in mm.
    pub slice_thickness: f64,
    /// Image orientation cosines (6 values: row_cos + col_cos).
    pub orientation: [f64; 6],
    /// Window center for display.
    pub window_center: f64,
    /// Window width for display.
    pub window_width: f64,
    /// Rescale slope (Hounsfield units: HU = pixel * slope + intercept).
    pub rescale_slope: f64,
    /// Rescale intercept.
    pub rescale_intercept: f64,
}
impl DicomImageData {
    /// Create DICOM image data.
    pub fn new(rows: usize, cols: usize, n_slices: usize) -> Self {
        Self {
            pixels: vec![0; rows * cols * n_slices],
            rows,
            cols,
            n_slices,
            pixel_spacing: [1.0, 1.0],
            slice_thickness: 1.0,
            orientation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            window_center: 40.0,
            window_width: 400.0,
            rescale_slope: 1.0,
            rescale_intercept: -1024.0,
        }
    }
    /// Convert pixel value to Hounsfield Units.
    pub fn to_hu(&self, pixel: u16) -> f64 {
        pixel as f64 * self.rescale_slope + self.rescale_intercept
    }
    /// Apply window/level mapping to HU value -> display [0, 255].
    pub fn window_level(&self, hu: f64) -> u8 {
        let low = self.window_center - self.window_width / 2.0;
        let high = self.window_center + self.window_width / 2.0;
        let norm = (hu - low) / (high - low);
        (norm.clamp(0.0, 1.0) * 255.0) as u8
    }
    /// Physical size in mm.
    pub fn physical_size(&self) -> [f64; 3] {
        [
            self.rows as f64 * self.pixel_spacing[0],
            self.cols as f64 * self.pixel_spacing[1],
            self.n_slices as f64 * self.slice_thickness,
        ]
    }
    /// Get voxel at (row, col, slice).
    pub fn voxel(&self, row: usize, col: usize, slice: usize) -> Option<u16> {
        if row < self.rows && col < self.cols && slice < self.n_slices {
            Some(self.pixels[slice * self.rows * self.cols + row * self.cols + col])
        } else {
            None
        }
    }
}
/// STL binary exporter.
pub struct StlExporter;
impl StlExporter {
    /// Serialise a list of triangles to a binary STL buffer (80-byte header + triangles).
    pub fn to_binary(triangles: &[StlTriangle]) -> Vec<u8> {
        let mut buf = vec![0u8; 80];
        let count = triangles.len() as u32;
        buf.extend_from_slice(&count.to_le_bytes());
        for tri in triangles {
            for &v in &tri.normal {
                buf.extend_from_slice(&v.to_le_bytes());
            }
            for vtx in &tri.vertices {
                for &c in vtx {
                    buf.extend_from_slice(&c.to_le_bytes());
                }
            }
            buf.extend_from_slice(&0u16.to_le_bytes());
        }
        buf
    }
    /// Parse a binary STL buffer back into a vector of triangles.
    pub fn from_binary(data: &[u8]) -> Vec<StlTriangle> {
        if data.len() < 84 {
            return Vec::new();
        }
        let count = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
        let mut triangles = Vec::with_capacity(count);
        let mut pos = 84usize;
        for _ in 0..count {
            if pos + 50 > data.len() {
                break;
            }
            let read_f32 = |p: &mut usize| -> f32 {
                let v = f32::from_le_bytes([data[*p], data[*p + 1], data[*p + 2], data[*p + 3]]);
                *p += 4;
                v
            };
            let nx = read_f32(&mut pos);
            let ny = read_f32(&mut pos);
            let nz = read_f32(&mut pos);
            let mut verts = [[0f32; 3]; 3];
            for vtx in &mut verts {
                vtx[0] = read_f32(&mut pos);
                vtx[1] = read_f32(&mut pos);
                vtx[2] = read_f32(&mut pos);
            }
            pos += 2;
            triangles.push(StlTriangle {
                normal: [nx, ny, nz],
                vertices: verts,
            });
        }
        triangles
    }
}
/// DICOM Value Representation (VR) types.
#[derive(Clone, Debug, PartialEq)]
pub enum DicomVr {
    /// Unsigned Short (2 bytes).
    Us(u16),
    /// Decimal String (ASCII numeric).
    Ds(f64),
    /// Long String (64-char max).
    Lo(String),
    /// Unique Identifier (UID string).
    Ui(String),
    /// Sequence of items (nested data sets).
    Sq(Vec<HashMap<DicomTag, DicomVr>>),
    /// Pixel Data (raw bytes).
    OW(Vec<u8>),
}
impl DicomVr {
    /// Return the VR type name as a 2-character DICOM code.
    pub fn vr_code(&self) -> &'static str {
        match self {
            DicomVr::Us(_) => "US",
            DicomVr::Ds(_) => "DS",
            DicomVr::Lo(_) => "LO",
            DicomVr::Ui(_) => "UI",
            DicomVr::Sq(_) => "SQ",
            DicomVr::OW(_) => "OW",
        }
    }
    /// Return the unsigned-short value if this is a US VR.
    pub fn as_us(&self) -> Option<u16> {
        if let DicomVr::Us(v) = self {
            Some(*v)
        } else {
            None
        }
    }
    /// Return the decimal string value if this is a DS VR.
    pub fn as_ds(&self) -> Option<f64> {
        if let DicomVr::Ds(v) = self {
            Some(*v)
        } else {
            None
        }
    }
    /// Return the string value for LO or UI VRs.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            DicomVr::Lo(s) | DicomVr::Ui(s) => Some(s.as_str()),
            _ => None,
        }
    }
}
/// MINC dimension descriptor (mirrors minc_volume dim structure).
#[derive(Clone, Debug)]
pub struct MincDimension {
    /// Dimension name (e.g. "xspace", "yspace", "zspace", "time").
    pub name: String,
    /// Number of voxels along this dimension.
    pub length: usize,
    /// Voxel step size in mm (can be negative for direction).
    pub step: f64,
    /// World-coordinate start value.
    pub start: f64,
}
impl MincDimension {
    /// Create a new MINC dimension.
    pub fn new(name: &str, length: usize, step: f64, start: f64) -> Self {
        Self {
            name: name.to_string(),
            length,
            step,
            start,
        }
    }
    /// Return the world coordinate of voxel index `i`.
    pub fn world_coord(&self, i: usize) -> f64 {
        self.start + i as f64 * self.step
    }
}
