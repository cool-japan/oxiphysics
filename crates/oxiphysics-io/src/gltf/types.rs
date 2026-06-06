//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::escape_json;

/// Light type for the `KHR_lights_punctual` extension.
#[derive(Debug, Clone, PartialEq)]
pub enum LightType {
    /// Omnidirectional point light.
    Point,
    /// Directional parallel light.
    Directional,
    /// Spot light with cone angles.
    Spot {
        /// Inner cone half-angle (radians).
        inner_cone_angle: f32,
        /// Outer cone half-angle (radians).
        outer_cone_angle: f32,
    },
}
/// Camera projection type.
#[derive(Debug, Clone)]
pub enum CameraType {
    /// Perspective projection.
    Perspective {
        /// Vertical field of view in radians.
        yfov: f32,
        /// Near clip plane.
        znear: f32,
        /// Far clip plane (None = infinite).
        zfar: Option<f32>,
        /// Aspect ratio (width / height). None = use viewport.
        aspect_ratio: Option<f32>,
    },
    /// Orthographic projection.
    Orthographic {
        /// Half-width of the view volume.
        xmag: f32,
        /// Half-height of the view volume.
        ymag: f32,
        /// Near clip plane.
        znear: f32,
        /// Far clip plane.
        zfar: f32,
    },
}
/// A validation issue found in a glTF scene.
#[derive(Debug, Clone)]
pub struct ValidationIssue {
    /// Human-readable description of the issue.
    pub message: String,
    /// True for errors (must fix), false for warnings (should fix).
    pub is_error: bool,
}
impl ValidationIssue {
    pub(crate) fn error(msg: impl Into<String>) -> Self {
        ValidationIssue {
            message: msg.into(),
            is_error: true,
        }
    }
    pub(crate) fn warning(msg: impl Into<String>) -> Self {
        ValidationIssue {
            message: msg.into(),
            is_error: false,
        }
    }
}
/// A single animation channel.
#[derive(Debug, Clone)]
pub struct GltfAnimationChannel {
    /// The channel target.
    pub target: GltfAnimationTarget,
    /// Interpolation method: "LINEAR", "STEP", "CUBICSPLINE".
    pub interpolation: String,
    /// Input (time) keyframe values.
    pub input_times: Vec<f32>,
    /// Output keyframe values (flattened; e.g., 3 floats per VEC3 keyframe).
    pub output_values: Vec<f32>,
}
/// A single keyframe channel for a joint.
#[derive(Debug, Clone)]
pub struct JointChannel {
    /// Index into the skeleton's joints array.
    pub joint_index: usize,
    /// Target path: "translation", "rotation", or "scale".
    pub path: String,
    /// Input (time) values.
    pub times: Vec<f32>,
    /// Output values (flattened; e.g., 3 floats/VEC3 for translation).
    pub values: Vec<f32>,
    /// Interpolation: "LINEAR", "STEP", or "CUBICSPLINE".
    pub interpolation: String,
}
/// Material properties for a glTF primitive.
#[derive(Debug, Clone)]
pub struct GltfMaterial {
    /// Material name.
    pub name: String,
    /// Base color factor \[r, g, b, a\] (PBR metallic-roughness).
    pub base_color_factor: [f32; 4],
    /// Metallic factor (0.0 = dielectric, 1.0 = metal).
    pub metallic_factor: f32,
    /// Roughness factor (0.0 = smooth, 1.0 = rough).
    pub roughness_factor: f32,
    /// Emissive factor \[r, g, b\].
    pub emissive_factor: [f32; 3],
    /// Alpha mode: "OPAQUE", "MASK", or "BLEND".
    pub alpha_mode: String,
    /// Alpha cutoff (only used when alpha_mode is "MASK").
    pub alpha_cutoff: f32,
    /// Whether the material is double-sided.
    pub double_sided: bool,
}
impl GltfMaterial {
    /// Check if the material is fully opaque (alpha = 1.0 and OPAQUE mode).
    pub fn is_opaque(&self) -> bool {
        self.alpha_mode == "OPAQUE" && (self.base_color_factor[3] - 1.0).abs() < 1e-6
    }
    /// Check if the material is emissive (any emissive factor > 0).
    pub fn is_emissive(&self) -> bool {
        self.emissive_factor.iter().any(|&v| v > 0.0)
    }
}
/// A joint (bone) in a skeleton hierarchy.
#[derive(Debug, Clone)]
pub struct Joint {
    /// Human-readable name.
    pub name: String,
    /// Translation \[x, y, z\].
    pub translation: [f64; 3],
    /// Rotation quaternion \[x, y, z, w\].
    pub rotation: [f64; 4],
    /// Scale \[x, y, z\].
    pub scale: [f64; 3],
    /// Child joint indices.
    pub children: Vec<usize>,
    /// Inverse bind matrix stored in column-major order (16 floats).
    pub inverse_bind_matrix: [f32; 16],
}
/// A light source (KHR_lights_punctual extension).
#[derive(Debug, Clone)]
pub struct SceneLight {
    /// Light name.
    pub name: String,
    /// Light type.
    pub light_type: LightType,
    /// Linear RGB color.
    pub color: [f32; 3],
    /// Intensity (candela for point/spot, lux for directional).
    pub intensity: f32,
}
impl SceneLight {
    /// Return the light type string as used in glTF JSON.
    pub fn type_string(&self) -> &str {
        match &self.light_type {
            LightType::Point => "point",
            LightType::Directional => "directional",
            LightType::Spot { .. } => "spot",
        }
    }
}
/// A single draw primitive: positions, normals, and indices.
pub struct GltfPrimitive {
    /// Per-vertex positions.
    pub positions: Vec<[f32; 3]>,
    /// Per-vertex normals.
    pub normals: Vec<[f32; 3]>,
    /// Per-vertex UV texture coordinates (TEXCOORD_0).
    pub texcoords: Vec<[f32; 2]>,
    /// Triangle indices.
    pub indices: Vec<u32>,
}
impl GltfPrimitive {
    /// Compute the axis-aligned bounding box of positions.
    /// Returns `(min, max)` as `[f32; 3]` arrays.
    pub fn bounding_box(&self) -> ([f32; 3], [f32; 3]) {
        let mut min = [f32::INFINITY; 3];
        let mut max = [f32::NEG_INFINITY; 3];
        for p in &self.positions {
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
    /// Return the number of triangles (indices / 3).
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
    /// Return the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }
    /// Extract triangle vertex positions as triplets.
    pub fn extract_triangles(&self) -> Vec<[[f32; 3]; 3]> {
        let mut tris = Vec::new();
        let mut i = 0;
        while i + 2 < self.indices.len() {
            let a = self.indices[i] as usize;
            let b = self.indices[i + 1] as usize;
            let c = self.indices[i + 2] as usize;
            if a < self.positions.len() && b < self.positions.len() && c < self.positions.len() {
                tris.push([self.positions[a], self.positions[b], self.positions[c]]);
            }
            i += 3;
        }
        tris
    }
}
/// A named mesh containing one or more primitives.
pub struct GltfMesh {
    /// Human-readable mesh name.
    pub name: String,
    /// List of draw primitives.
    pub primitives: Vec<GltfPrimitive>,
}
impl GltfMesh {
    /// Total vertex count across all primitives.
    pub fn total_vertex_count(&self) -> usize {
        self.primitives.iter().map(|p| p.vertex_count()).sum()
    }
    /// Total triangle count across all primitives.
    pub fn total_triangle_count(&self) -> usize {
        self.primitives.iter().map(|p| p.triangle_count()).sum()
    }
}
/// Element type (SCALAR, VEC2, VEC3, VEC4, MAT4, …).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AccessorType {
    /// Single scalar value.
    Scalar,
    /// 2-component vector.
    Vec2,
    /// 3-component vector.
    Vec3,
    /// 4-component vector.
    Vec4,
    /// 2×2 column-major matrix.
    Mat2,
    /// 3×3 column-major matrix.
    Mat3,
    /// 4×4 column-major matrix.
    Mat4,
}
impl AccessorType {
    /// Number of components per element.
    pub fn num_components(self) -> usize {
        match self {
            AccessorType::Scalar => 1,
            AccessorType::Vec2 => 2,
            AccessorType::Vec3 => 3,
            AccessorType::Vec4 => 4,
            AccessorType::Mat2 => 4,
            AccessorType::Mat3 => 9,
            AccessorType::Mat4 => 16,
        }
    }
    /// glTF JSON string representation.
    pub fn as_str(self) -> &'static str {
        match self {
            AccessorType::Scalar => "SCALAR",
            AccessorType::Vec2 => "VEC2",
            AccessorType::Vec3 => "VEC3",
            AccessorType::Vec4 => "VEC4",
            AccessorType::Mat2 => "MAT2",
            AccessorType::Mat3 => "MAT3",
            AccessorType::Mat4 => "MAT4",
        }
    }
}
/// A skeleton animation: name + per-joint channels.
#[derive(Debug, Clone)]
pub struct SkeletonAnimation {
    /// Animation name.
    pub name: String,
    /// Per-joint keyframe channels.
    pub joint_channels: Vec<JointChannel>,
}
impl SkeletonAnimation {
    /// Duration of the animation (max time across all channels).
    pub fn duration(&self) -> f32 {
        self.joint_channels
            .iter()
            .filter_map(|c| c.times.last().copied())
            .fold(0.0f32, f32::max)
    }
    /// Total keyframe count across all channels.
    pub fn total_keyframes(&self) -> usize {
        self.joint_channels.iter().map(|c| c.times.len()).sum()
    }
}
/// A buffer view defines a slice of a binary buffer.
#[derive(Debug, Clone)]
pub struct GltfBufferView {
    /// Index into the buffers array.
    pub buffer: usize,
    /// Byte offset into the buffer.
    pub byte_offset: usize,
    /// Length in bytes.
    pub byte_length: usize,
    /// Optional byte stride (for interleaved data).
    pub byte_stride: Option<usize>,
}
/// A typed array view over a raw byte buffer.
///
/// Mirrors the glTF `accessor` concept — a window into a `bufferView` with
/// a specific element type and component type.
pub struct TypedAccessor {
    /// Human-readable label.
    pub name: String,
    /// Index of the owning buffer view.
    pub buffer_view_index: usize,
    /// Byte offset within the buffer view.
    pub byte_offset: usize,
    /// Component type (FLOAT, UNSIGNED_INT, …).
    pub component_type: ComponentType,
    /// Element type (SCALAR, VEC3, …).
    pub accessor_type: AccessorType,
    /// Number of elements (not bytes, not components).
    pub count: usize,
    /// Optional min values per component (length = num_components).
    pub min_values: Vec<f64>,
    /// Optional max values per component.
    pub max_values: Vec<f64>,
}
impl TypedAccessor {
    /// Create a new accessor.
    pub fn new(
        name: impl Into<String>,
        buffer_view_index: usize,
        byte_offset: usize,
        component_type: ComponentType,
        accessor_type: AccessorType,
        count: usize,
    ) -> Self {
        Self {
            name: name.into(),
            buffer_view_index,
            byte_offset,
            component_type,
            accessor_type,
            count,
            min_values: Vec::new(),
            max_values: Vec::new(),
        }
    }
    /// Total byte length of this accessor's data.
    pub fn byte_length(&self) -> usize {
        self.count * self.accessor_type.num_components() * self.component_type.byte_size()
    }
    /// Serialize to a glTF JSON object fragment (returns a JSON object string).
    pub fn to_json(&self) -> String {
        format!(
            r#"{{ "bufferView": {bv}, "byteOffset": {bo}, "componentType": {ct}, "type": "{at}", "count": {c} }}"#,
            bv = self.buffer_view_index,
            bo = self.byte_offset,
            ct = self.component_type.component_type_code(),
            at = self.accessor_type.as_str(),
            c = self.count,
        )
    }
    /// Decode this accessor's data as Vec`f32` from a raw byte buffer.
    ///
    /// Only valid when `component_type == Float`.
    pub fn decode_f32(&self, buffer: &[u8]) -> Option<Vec<f32>> {
        if self.component_type != ComponentType::Float {
            return None;
        }
        let start = self.byte_offset;
        let len = self.byte_length();
        let end = start + len;
        if end > buffer.len() {
            return None;
        }
        let slice = &buffer[start..end];
        let floats: Vec<f32> = slice
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect();
        Some(floats)
    }
    /// Decode as Vec`u32` indices (only valid for UNSIGNED_INT accessors).
    pub fn decode_u32(&self, buffer: &[u8]) -> Option<Vec<u32>> {
        if self.component_type != ComponentType::UnsignedInt {
            return None;
        }
        let start = self.byte_offset;
        let len = self.byte_length();
        let end = start + len;
        if end > buffer.len() {
            return None;
        }
        let slice = &buffer[start..end];
        let indices: Vec<u32> = slice
            .chunks_exact(4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            .collect();
        Some(indices)
    }
}
/// An animation channel target.
#[derive(Debug, Clone)]
pub struct GltfAnimationTarget {
    /// Index of the target node.
    pub node: usize,
    /// Path: "translation", "rotation", "scale", or "weights".
    pub path: String,
}
/// A node in a glTF scene graph.
pub struct GltfNode {
    /// Human-readable name of the node.
    pub name: String,
    /// Optional index into the meshes array.
    pub mesh: Option<usize>,
    /// Translation vector \[x, y, z\].
    pub translation: [f64; 3],
    /// Rotation quaternion \[x, y, z, w\].
    pub rotation: [f64; 4],
    /// Scale vector \[x, y, z\].
    pub scale: [f64; 3],
    /// Child node indices for scene graph traversal.
    pub children: Vec<usize>,
}
/// Interpolation mode for animation channels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Interpolation {
    /// Linear interpolation between keyframes.
    Linear,
    /// Hold value until next keyframe.
    Step,
    /// Cubic spline (requires tangents).
    CubicSpline,
}
impl Interpolation {
    /// glTF JSON string.
    pub fn as_str(&self) -> &'static str {
        match self {
            Interpolation::Linear => "LINEAR",
            Interpolation::Step => "STEP",
            Interpolation::CubicSpline => "CUBICSPLINE",
        }
    }
}
/// A skeleton: a collection of joints.
pub struct Skeleton {
    /// List of joints.
    pub joints: Vec<Joint>,
}
impl Skeleton {
    /// Create an empty skeleton.
    pub fn new() -> Self {
        Skeleton { joints: Vec::new() }
    }
    /// Add a joint and return its index.
    pub fn add_joint(&mut self, joint: Joint) -> usize {
        let idx = self.joints.len();
        self.joints.push(joint);
        idx
    }
    /// Return the number of joints.
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }
    /// Export a `SkeletonAnimation` as a simplified glTF animation JSON object.
    ///
    /// Produces a JSON string suitable for embedding in a `"animations"` array.
    pub fn export_animation_json(&self, anim: &SkeletonAnimation) -> String {
        let mut out = String::new();
        out.push_str("{\n");
        out.push_str(&format!("  \"name\": \"{}\",\n", escape_json(&anim.name)));
        out.push_str("  \"channels\": [\n");
        for (ci, ch) in anim.joint_channels.iter().enumerate() {
            let joint_name = self
                .joints
                .get(ch.joint_index)
                .map(|j| j.name.as_str())
                .unwrap_or("unknown");
            out.push_str("    {\n");
            out.push_str(&format!(
                "      \"joint\": \"{}\",\n",
                escape_json(joint_name)
            ));
            out.push_str(&format!("      \"path\": \"{}\",\n", escape_json(&ch.path)));
            out.push_str(&format!(
                "      \"interpolation\": \"{}\",\n",
                escape_json(&ch.interpolation)
            ));
            let times_str: Vec<String> = ch.times.iter().map(|t| format!("{t}")).collect();
            out.push_str(&format!("      \"times\": [{}]\n", times_str.join(", ")));
            if ci + 1 < anim.joint_channels.len() {
                out.push_str("    },\n");
            } else {
                out.push_str("    }\n");
            }
        }
        out.push_str("  ]\n");
        out.push('}');
        out
    }
}
/// A camera in the scene.
#[derive(Debug, Clone)]
pub struct SceneCamera {
    /// Camera name.
    pub name: String,
    /// Projection type.
    pub camera_type: CameraType,
}
/// Component type codes used in glTF accessors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentType {
    /// Unsigned byte (1 byte).
    UnsignedByte = 5121,
    /// Unsigned short (2 bytes).
    UnsignedShort = 5123,
    /// Unsigned int (4 bytes).
    UnsignedInt = 5125,
    /// Float (4 bytes).
    Float = 5126,
}
impl ComponentType {
    /// Byte size of one component.
    pub fn byte_size(self) -> usize {
        match self {
            ComponentType::UnsignedByte => 1,
            ComponentType::UnsignedShort => 2,
            ComponentType::UnsignedInt => 4,
            ComponentType::Float => 4,
        }
    }
    /// glTF JSON numeric value.
    pub fn component_type_code(self) -> u32 {
        self as u32
    }
}
/// Builder for a PBR metallic-roughness material description.
///
/// Generates a glTF-compatible JSON material object fragment.
#[derive(Debug, Clone)]
pub struct PbrMaterialBuilder {
    /// Material name.
    pub name: String,
    /// Base colour factor \[R, G, B, A\] in linear sRGB.
    pub base_color_factor: [f32; 4],
    /// Metallic factor in \[0.0, 1.0\].
    pub metallic_factor: f32,
    /// Roughness factor in \[0.0, 1.0\].
    pub roughness_factor: f32,
    /// Emissive colour factor \[R, G, B\].
    pub emissive_factor: [f32; 3],
    /// Double-sided rendering.
    pub double_sided: bool,
    /// Alpha mode: "OPAQUE", "MASK", or "BLEND".
    pub alpha_mode: String,
    /// Alpha cutoff (for MASK mode).
    pub alpha_cutoff: f32,
}
impl PbrMaterialBuilder {
    /// Create a new builder with default values.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }
    /// Set the base colour factor.
    pub fn base_color(mut self, r: f32, g: f32, b: f32, a: f32) -> Self {
        self.base_color_factor = [r, g, b, a];
        self
    }
    /// Set metallic and roughness factors.
    pub fn metallic_roughness(mut self, metallic: f32, roughness: f32) -> Self {
        self.metallic_factor = metallic;
        self.roughness_factor = roughness;
        self
    }
    /// Set emissive factor.
    pub fn emissive(mut self, r: f32, g: f32, b: f32) -> Self {
        self.emissive_factor = [r, g, b];
        self
    }
    /// Enable double-sided rendering.
    pub fn double_sided(mut self, ds: bool) -> Self {
        self.double_sided = ds;
        self
    }
    /// Set alpha mode.
    pub fn alpha_mode(mut self, mode: impl Into<String>) -> Self {
        self.alpha_mode = mode.into();
        self
    }
    /// Build as a JSON string (glTF material object).
    pub fn to_json(&self) -> String {
        let bc = &self.base_color_factor;
        let em = &self.emissive_factor;
        format!(
            r#"{{
  "name": "{name}",
  "pbrMetallicRoughness": {{
    "baseColorFactor": [{r:.6}, {g:.6}, {b:.6}, {a:.6}],
    "metallicFactor": {mf:.6},
    "roughnessFactor": {rf:.6}
  }},
  "emissiveFactor": [{er:.6}, {eg:.6}, {eb:.6}],
  "doubleSided": {ds},
  "alphaMode": "{am}",
  "alphaCutoff": {ac:.6}
}}"#,
            name = self.name,
            r = bc[0],
            g = bc[1],
            b = bc[2],
            a = bc[3],
            mf = self.metallic_factor,
            rf = self.roughness_factor,
            er = em[0],
            eg = em[1],
            eb = em[2],
            ds = self.double_sided,
            am = self.alpha_mode,
            ac = self.alpha_cutoff,
        )
    }
    /// Convert to a [`GltfMaterial`] instance (for use in a [`GltfScene`]).
    pub fn build(&self) -> GltfMaterial {
        GltfMaterial {
            name: self.name.clone(),
            base_color_factor: self.base_color_factor,
            metallic_factor: self.metallic_factor,
            roughness_factor: self.roughness_factor,
            emissive_factor: self.emissive_factor,
            double_sided: self.double_sided,
            alpha_mode: self.alpha_mode.clone(),
            alpha_cutoff: self.alpha_cutoff,
        }
    }
}
/// Builder for a single animation channel (translation / rotation / scale).
pub struct AnimationChannelBuilder {
    /// Target node index.
    pub node_index: usize,
    /// Target path: "translation", "rotation", or "scale".
    pub path: String,
    /// Interpolation mode.
    pub interpolation: Interpolation,
    /// Ordered keyframes (should be sorted by time).
    pub keyframes: Vec<Keyframe>,
}
impl AnimationChannelBuilder {
    /// Create a new channel builder.
    pub fn new(node_index: usize, path: impl Into<String>, interpolation: Interpolation) -> Self {
        Self {
            node_index,
            path: path.into(),
            interpolation,
            keyframes: Vec::new(),
        }
    }
    /// Append a keyframe.
    pub fn push(mut self, time: f32, value: Vec<f32>) -> Self {
        self.keyframes.push(Keyframe::new(time, value));
        self
    }
    /// Duration of the channel (time of the last keyframe, or 0 if empty).
    pub fn duration(&self) -> f32 {
        self.keyframes.last().map(|k| k.time).unwrap_or(0.0)
    }
    /// Number of keyframes.
    pub fn len(&self) -> usize {
        self.keyframes.len()
    }
    /// True if there are no keyframes.
    pub fn is_empty(&self) -> bool {
        self.keyframes.is_empty()
    }
    /// Collect all time stamps.
    pub fn times(&self) -> Vec<f32> {
        self.keyframes.iter().map(|k| k.time).collect()
    }
    /// Collect all value payloads flattened into one Vec.
    pub fn values_flat(&self) -> Vec<f32> {
        self.keyframes
            .iter()
            .flat_map(|k| k.value.iter().copied())
            .collect()
    }
    /// Serialize sampler + channel pair to JSON fragments.
    pub fn to_json_fragments(
        &self,
        sampler_index: usize,
        times_accessor: usize,
        values_accessor: usize,
    ) -> (String, String) {
        let sampler = format!(
            r#"{{ "input": {ti}, "interpolation": "{interp}", "output": {vi} }}"#,
            ti = times_accessor,
            interp = self.interpolation.as_str(),
            vi = values_accessor,
        );
        let channel = format!(
            r#"{{ "sampler": {si}, "target": {{ "node": {ni}, "path": "{path}" }} }}"#,
            si = sampler_index,
            ni = self.node_index,
            path = self.path,
        );
        (sampler, channel)
    }
}
/// A named animation containing multiple channels.
#[derive(Debug, Clone)]
pub struct GltfAnimation {
    /// Animation name.
    pub name: String,
    /// Animation channels.
    pub channels: Vec<GltfAnimationChannel>,
}
impl GltfAnimation {
    /// Total duration of the animation (max input time across all channels).
    pub fn duration(&self) -> f32 {
        self.channels
            .iter()
            .filter_map(|ch| ch.input_times.last().copied())
            .fold(0.0f32, f32::max)
    }
    /// Total number of keyframes across all channels.
    pub fn total_keyframe_count(&self) -> usize {
        self.channels.iter().map(|ch| ch.input_times.len()).sum()
    }
}
/// A single keyframe with a time stamp and a value payload.
#[derive(Debug, Clone)]
pub struct Keyframe {
    /// Time in seconds.
    pub time: f32,
    /// Value payload (length depends on the target path).
    pub value: Vec<f32>,
}
impl Keyframe {
    /// Create a new keyframe.
    pub fn new(time: f32, value: Vec<f32>) -> Self {
        Self { time, value }
    }
}
/// An accessor describes how to read typed data from a buffer view.
#[derive(Debug, Clone)]
pub struct GltfAccessor {
    /// Index into the buffer views array.
    pub buffer_view: usize,
    /// Byte offset within the buffer view.
    pub byte_offset: usize,
    /// Component type (5126 = FLOAT, 5123 = UNSIGNED_SHORT, 5125 = UNSIGNED_INT).
    pub component_type: u32,
    /// Number of elements.
    pub count: usize,
    /// Type string ("SCALAR", "VEC2", "VEC3", "VEC4", "MAT4").
    pub element_type: String,
}
impl GltfAccessor {
    /// Number of components per element.
    pub fn components_per_element(&self) -> usize {
        match self.element_type.as_str() {
            "SCALAR" => 1,
            "VEC2" => 2,
            "VEC3" => 3,
            "VEC4" => 4,
            "MAT2" => 4,
            "MAT3" => 9,
            "MAT4" => 16,
            _ => 1,
        }
    }
    /// Size of a single component in bytes.
    pub fn component_size(&self) -> usize {
        match self.component_type {
            5120 => 1,
            5121 => 1,
            5122 => 2,
            5123 => 2,
            5125 => 4,
            5126 => 4,
            _ => 4,
        }
    }
    /// Total byte length of data described by this accessor.
    pub fn byte_length(&self) -> usize {
        self.count * self.components_per_element() * self.component_size()
    }
}
/// A morph target (blend shape) that stores position deltas.
#[derive(Debug, Clone)]
pub struct MorphTarget {
    /// Human-readable name of the morph target.
    pub name: String,
    /// Per-vertex position deltas applied at weight = 1.0.
    pub position_deltas: Vec<[f32; 3]>,
}
impl MorphTarget {
    /// Apply this morph target to a set of base positions with the given weight.
    ///
    /// Returns a new `Vec` with blended positions.
    pub fn apply(&self, base_positions: &[[f32; 3]], weight: f32) -> Vec<[f32; 3]> {
        let n = base_positions.len().min(self.position_deltas.len());
        let mut out = base_positions.to_vec();
        for (dst, delta) in out[..n].iter_mut().zip(self.position_deltas[..n].iter()) {
            dst[0] += delta[0] * weight;
            dst[1] += delta[1] * weight;
            dst[2] += delta[2] * weight;
        }
        out
    }
}
/// A glTF 2.0 scene: a collection of nodes and meshes.
pub struct GltfScene {
    /// Scene nodes.
    pub nodes: Vec<GltfNode>,
    /// Scene meshes.
    pub meshes: Vec<GltfMesh>,
    /// Accessors.
    pub accessors: Vec<GltfAccessor>,
    /// Buffer views.
    pub buffer_views: Vec<GltfBufferView>,
    /// Animations.
    pub animations: Vec<GltfAnimation>,
    /// Materials.
    pub materials: Vec<GltfMaterial>,
    /// Cameras.
    pub cameras: Vec<SceneCamera>,
    /// Lights (KHR_lights_punctual).
    pub lights: Vec<SceneLight>,
}
impl GltfScene {
    /// Create an empty scene.
    pub fn new() -> Self {
        GltfScene {
            nodes: Vec::new(),
            meshes: Vec::new(),
            accessors: Vec::new(),
            buffer_views: Vec::new(),
            animations: Vec::new(),
            materials: Vec::new(),
            cameras: Vec::new(),
            lights: Vec::new(),
        }
    }
    /// Add a camera and return its index.
    pub fn add_camera(&mut self, cam: SceneCamera) -> usize {
        let idx = self.cameras.len();
        self.cameras.push(cam);
        idx
    }
    /// Add a light and return its index.
    pub fn add_light(&mut self, light: SceneLight) -> usize {
        let idx = self.lights.len();
        self.lights.push(light);
        idx
    }
    /// Serialize the scene as a glTF 2.0 JSON string including cameras and lights.
    ///
    /// Cameras are written in the `cameras` array.
    /// Lights use the `KHR_lights_punctual` extension.
    pub fn to_json_with_hierarchy(&self) -> String {
        let mut out = self.to_json();
        if out.ends_with('}') {
            out.pop();
        }
        out.push_str(",\n");
        if !self.cameras.is_empty() {
            out.push_str("  \"cameras\": [\n");
            for (ci, cam) in self.cameras.iter().enumerate() {
                out.push_str("    {\n");
                out.push_str(&format!(
                    "      \"name\": \"{}\",\n",
                    escape_json(&cam.name)
                ));
                match &cam.camera_type {
                    CameraType::Perspective {
                        yfov,
                        znear,
                        zfar,
                        aspect_ratio,
                    } => {
                        out.push_str("      \"type\": \"perspective\",\n");
                        out.push_str("      \"perspective\": {\n");
                        out.push_str(&format!("        \"yfov\": {},\n", yfov));
                        out.push_str(&format!("        \"znear\": {}", znear));
                        if let Some(zf) = zfar {
                            out.push_str(&format!(",\n        \"zfar\": {}", zf));
                        }
                        if let Some(ar) = aspect_ratio {
                            out.push_str(&format!(",\n        \"aspectRatio\": {}", ar));
                        }
                        out.push_str("\n      }\n");
                    }
                    CameraType::Orthographic {
                        xmag,
                        ymag,
                        znear,
                        zfar,
                    } => {
                        out.push_str("      \"type\": \"orthographic\",\n");
                        out.push_str("      \"orthographic\": {\n");
                        out.push_str(&format!("        \"xmag\": {},\n", xmag));
                        out.push_str(&format!("        \"ymag\": {},\n", ymag));
                        out.push_str(&format!("        \"znear\": {},\n", znear));
                        out.push_str(&format!("        \"zfar\": {}\n", zfar));
                        out.push_str("      }\n");
                    }
                }
                if ci + 1 < self.cameras.len() {
                    out.push_str("    },\n");
                } else {
                    out.push_str("    }\n");
                }
            }
            out.push_str("  ],\n");
        }
        if !self.lights.is_empty() {
            out.push_str("  \"extensions\": {\n");
            out.push_str("    \"KHR_lights_punctual\": {\n");
            out.push_str("      \"lights\": [\n");
            for (li, light) in self.lights.iter().enumerate() {
                out.push_str("        {\n");
                out.push_str(&format!(
                    "          \"name\": \"{}\",\n",
                    escape_json(&light.name)
                ));
                out.push_str(&format!(
                    "          \"type\": \"{}\",\n",
                    light.type_string()
                ));
                out.push_str(&format!(
                    "          \"color\": [{}, {}, {}],\n",
                    light.color[0], light.color[1], light.color[2]
                ));
                out.push_str(&format!("          \"intensity\": {}", light.intensity));
                if let LightType::Spot {
                    inner_cone_angle,
                    outer_cone_angle,
                } = &light.light_type
                {
                    out.push_str(
                        &format!(
                            ",\n          \"spot\": {{ \"innerConeAngle\": {}, \"outerConeAngle\": {} }}",
                            inner_cone_angle, outer_cone_angle
                        ),
                    );
                }
                out.push('\n');
                if li + 1 < self.lights.len() {
                    out.push_str("        },\n");
                } else {
                    out.push_str("        }\n");
                }
            }
            out.push_str("      ]\n");
            out.push_str("    }\n");
            out.push_str("  }\n");
        } else if out.ends_with(",\n") {
            out.truncate(out.len() - 2);
            out.push('\n');
        }
        out.push('}');
        out
    }
    /// Add a mesh to the scene and return its index.
    pub fn add_mesh(&mut self, mesh: GltfMesh) -> usize {
        let idx = self.meshes.len();
        self.meshes.push(mesh);
        idx
    }
    /// Add a node to the scene.
    pub fn add_node(&mut self, node: GltfNode) {
        self.nodes.push(node);
    }
    /// Return the number of meshes in the scene.
    pub fn mesh_count(&self) -> usize {
        self.meshes.len()
    }
    /// Return the number of nodes in the scene.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Traverse the scene graph depth-first, calling `visitor` with
    /// `(node_index, depth, accumulated_translation)` for each node.
    pub fn traverse_depth_first<F>(&self, mut visitor: F)
    where
        F: FnMut(usize, usize, [f64; 3]),
    {
        let mut is_child = vec![false; self.nodes.len()];
        for node in &self.nodes {
            for &child_idx in &node.children {
                if child_idx < is_child.len() {
                    is_child[child_idx] = true;
                }
            }
        }
        let mut stack: Vec<(usize, usize, [f64; 3])> = Vec::new();
        for i in (0..self.nodes.len()).rev() {
            if !is_child[i] {
                stack.push((i, 0, [0.0, 0.0, 0.0]));
            }
        }
        while let Some((idx, depth, parent_translation)) = stack.pop() {
            if idx >= self.nodes.len() {
                continue;
            }
            let node = &self.nodes[idx];
            let accumulated = [
                parent_translation[0] + node.translation[0],
                parent_translation[1] + node.translation[1],
                parent_translation[2] + node.translation[2],
            ];
            visitor(idx, depth, accumulated);
            for &child_idx in node.children.iter().rev() {
                stack.push((child_idx, depth + 1, accumulated));
            }
        }
    }
    /// Collect all mesh primitives in the scene with their associated node names.
    pub fn collect_mesh_primitives(&self) -> Vec<(&str, &GltfPrimitive)> {
        let mut result = Vec::new();
        for node in &self.nodes {
            if let Some(mesh_idx) = node.mesh
                && let Some(mesh) = self.meshes.get(mesh_idx)
            {
                for prim in &mesh.primitives {
                    result.push((node.name.as_str(), prim));
                }
            }
        }
        result
    }
    /// Add a material and return its index.
    pub fn add_material(&mut self, mat: GltfMaterial) -> usize {
        let idx = self.materials.len();
        self.materials.push(mat);
        idx
    }
    /// Add an animation and return its index.
    pub fn add_animation(&mut self, anim: GltfAnimation) -> usize {
        let idx = self.animations.len();
        self.animations.push(anim);
        idx
    }
    /// Find all nodes that reference a given mesh index.
    pub fn nodes_using_mesh(&self, mesh_idx: usize) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                if n.mesh == Some(mesh_idx) {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Return the total vertex count across all meshes.
    pub fn total_vertex_count(&self) -> usize {
        self.meshes.iter().map(|m| m.total_vertex_count()).sum()
    }
    /// Return the total triangle count across all meshes.
    pub fn total_triangle_count(&self) -> usize {
        self.meshes.iter().map(|m| m.total_triangle_count()).sum()
    }
    /// Serialize the scene as a glTF 2.0 JSON string.
    ///
    /// Accessor and bufferView entries are written as stubs (no binary buffer).
    pub fn to_json(&self) -> String {
        let mut out = String::new();
        out.push_str("{\n");
        out.push_str("  \"asset\": {\n");
        out.push_str("    \"version\": \"2.0\",\n");
        out.push_str("    \"generator\": \"OxiPhysics glTF writer\"\n");
        out.push_str("  },\n");
        out.push_str("  \"scene\": 0,\n");
        let node_indices: Vec<String> = (0..self.nodes.len()).map(|i| i.to_string()).collect();
        out.push_str("  \"scenes\": [\n    {\n      \"nodes\": [");
        out.push_str(&node_indices.join(", "));
        out.push_str("]\n    }\n  ],\n");
        out.push_str("  \"nodes\": [\n");
        for (i, node) in self.nodes.iter().enumerate() {
            out.push_str("    {\n");
            out.push_str(&format!(
                "      \"name\": \"{}\",\n",
                escape_json(&node.name)
            ));
            if let Some(mesh_idx) = node.mesh {
                out.push_str(&format!("      \"mesh\": {},\n", mesh_idx));
            }
            if !node.children.is_empty() {
                let children_str: Vec<String> =
                    node.children.iter().map(|c| c.to_string()).collect();
                out.push_str(&format!(
                    "      \"children\": [{}],\n",
                    children_str.join(", ")
                ));
            }
            out.push_str(&format!(
                "      \"translation\": [{}, {}, {}],\n",
                node.translation[0], node.translation[1], node.translation[2]
            ));
            out.push_str(&format!(
                "      \"rotation\": [{}, {}, {}, {}],\n",
                node.rotation[0], node.rotation[1], node.rotation[2], node.rotation[3]
            ));
            out.push_str(&format!(
                "      \"scale\": [{}, {}, {}]\n",
                node.scale[0], node.scale[1], node.scale[2]
            ));
            if i + 1 < self.nodes.len() {
                out.push_str("    },\n");
            } else {
                out.push_str("    }\n");
            }
        }
        out.push_str("  ],\n");
        out.push_str("  \"meshes\": [\n");
        for (mi, mesh) in self.meshes.iter().enumerate() {
            out.push_str("    {\n");
            out.push_str(&format!(
                "      \"name\": \"{}\",\n",
                escape_json(&mesh.name)
            ));
            out.push_str("      \"primitives\": [\n");
            for (pi, _prim) in mesh.primitives.iter().enumerate() {
                let pos_acc = pi * 3;
                let nrm_acc = pi * 3 + 1;
                let idx_acc = pi * 3 + 2;
                out.push_str("        {\n");
                out.push_str("          \"attributes\": {\n");
                out.push_str(&format!("            \"POSITION\": {},\n", pos_acc));
                out.push_str(&format!("            \"NORMAL\": {}\n", nrm_acc));
                out.push_str("          },\n");
                out.push_str(&format!("          \"indices\": {}\n", idx_acc));
                if pi + 1 < mesh.primitives.len() {
                    out.push_str("        },\n");
                } else {
                    out.push_str("        }\n");
                }
            }
            out.push_str("      ]\n");
            if mi + 1 < self.meshes.len() {
                out.push_str("    },\n");
            } else {
                out.push_str("    }\n");
            }
        }
        out.push_str("  ],\n");
        if !self.materials.is_empty() {
            out.push_str("  \"materials\": [\n");
            for (mi, mat) in self.materials.iter().enumerate() {
                out.push_str("    {\n");
                out.push_str(&format!(
                    "      \"name\": \"{}\",\n",
                    escape_json(&mat.name)
                ));
                out.push_str("      \"pbrMetallicRoughness\": {\n");
                out.push_str(&format!(
                    "        \"baseColorFactor\": [{}, {}, {}, {}],\n",
                    mat.base_color_factor[0],
                    mat.base_color_factor[1],
                    mat.base_color_factor[2],
                    mat.base_color_factor[3]
                ));
                out.push_str(&format!(
                    "        \"metallicFactor\": {},\n",
                    mat.metallic_factor
                ));
                out.push_str(&format!(
                    "        \"roughnessFactor\": {}\n",
                    mat.roughness_factor
                ));
                out.push_str("      },\n");
                out.push_str(&format!(
                    "      \"emissiveFactor\": [{}, {}, {}],\n",
                    mat.emissive_factor[0], mat.emissive_factor[1], mat.emissive_factor[2]
                ));
                out.push_str(&format!("      \"alphaMode\": \"{}\",\n", mat.alpha_mode));
                out.push_str(&format!("      \"doubleSided\": {}\n", mat.double_sided));
                if mi + 1 < self.materials.len() {
                    out.push_str("    },\n");
                } else {
                    out.push_str("    }\n");
                }
            }
            out.push_str("  ],\n");
        }
        out.push_str("  \"accessors\": [],\n");
        out.push_str("  \"bufferViews\": [],\n");
        out.push_str("  \"buffers\": []\n");
        out.push('}');
        out
    }
}
/// A mesh primitive extended with morph targets and blend weights.
pub struct MorphPrimitive {
    /// The base geometry.
    pub base: GltfPrimitive,
    /// List of morph targets.
    pub targets: Vec<MorphTarget>,
    /// Current blend weights (one per target).
    pub weights: Vec<f32>,
}
impl MorphPrimitive {
    /// Set the blend weight for a given target index (clamped to \[0, 1\]).
    pub fn set_weight(&mut self, target_idx: usize, weight: f32) {
        if target_idx < self.weights.len() {
            self.weights[target_idx] = weight.clamp(0.0, 1.0);
        }
    }
    /// Compute the blended positions using all target weights.
    pub fn blend(&self) -> Vec<[f32; 3]> {
        let mut result = self.base.positions.clone();
        for (target, &w) in self.targets.iter().zip(self.weights.iter()) {
            let n = result.len().min(target.position_deltas.len());
            for (dst, delta) in result[..n]
                .iter_mut()
                .zip(target.position_deltas[..n].iter())
            {
                dst[0] += delta[0] * w;
                dst[1] += delta[1] * w;
                dst[2] += delta[2] * w;
            }
        }
        result
    }
}
/// Writer that serialises a `GltfScene` into the GLB binary format.
///
/// GLB layout:
/// - 12-byte header: magic (4), version (4), total_length (4)
/// - JSON chunk:     chunk_length (4), chunk_type "JSON" (4), JSON payload (aligned to 4 bytes)
/// - BIN chunk (if any binary data): chunk_length (4), chunk_type "BIN\0" (4), data
pub struct GlbWriter {
    /// Whether to include a (empty) BIN chunk even when there is no binary data.
    pub include_empty_bin: bool,
}
impl GlbWriter {
    /// Create a new `GlbWriter`.
    pub fn new() -> Self {
        GlbWriter {
            include_empty_bin: false,
        }
    }

    /// Serialise the scene into GLB bytes (JSON-only, no vertex binary data).
    ///
    /// For a full binary GLB with packed vertex data use [`write_glb`](Self::write_glb).
    pub fn write(&self, scene: &GltfScene) -> Vec<u8> {
        let json = scene.to_json();
        let json_bytes = json.as_bytes();
        let json_len = json_bytes.len();
        let json_padded_len = (json_len + 3) & !3;
        let json_padding = json_padded_len - json_len;
        let chunk_header_size = 8usize;
        let header_size = 12usize;
        let total_len = header_size
            + chunk_header_size
            + json_padded_len
            + if self.include_empty_bin {
                chunk_header_size
            } else {
                0
            };
        let mut out = Vec::with_capacity(total_len);
        out.extend_from_slice(b"glTF");
        out.extend_from_slice(&2u32.to_le_bytes());
        out.extend_from_slice(&(total_len as u32).to_le_bytes());
        out.extend_from_slice(&(json_padded_len as u32).to_le_bytes());
        out.extend_from_slice(&0x4E4F534Au32.to_le_bytes());
        out.extend_from_slice(json_bytes);
        out.extend(std::iter::repeat_n(0x20u8, json_padding));
        if self.include_empty_bin {
            out.extend_from_slice(&0u32.to_le_bytes());
            out.extend_from_slice(&0x004E4942u32.to_le_bytes());
        }
        out
    }

    /// Serialise the first mesh primitive in the scene as a proper glTF 2.0 Binary (GLB) file.
    ///
    /// GLB binary layout:
    /// 1. **12-byte header**: magic `glTF`, version=2, total_length
    /// 2. **JSON chunk**: chunk_length (u32 LE), chunk_type=`JSON` (0x4E4F534A), padded UTF-8 JSON
    /// 3. **BIN chunk**: chunk_length (u32 LE), chunk_type=`BIN\0` (0x004E4942), vertex/index data
    ///
    /// Binary buffer layout (flat, tightly packed, index buffer 4-byte aligned):
    /// - `positions` : `[f32; 3]` × n_vertices
    /// - `normals`   : `[f32; 3]` × n_vertices
    /// - `texcoords` : `[f32; 2]` × n_vertices (zero-filled if `primitive.texcoords` is empty)
    /// - `indices`   : `u32` × n_indices (offset aligned to 4 bytes)
    ///
    /// Returns an empty `Vec<u8>` if the scene has no meshes.
    pub fn write_glb(&self, scene: &GltfScene) -> Vec<u8> {
        // Collect all primitives across all meshes.
        let primitives: Vec<&GltfPrimitive> = scene
            .meshes
            .iter()
            .flat_map(|m| m.primitives.iter())
            .collect();

        if primitives.is_empty() {
            return Vec::new();
        }

        // Build flat binary buffer: positions | normals | texcoords | (align) | indices
        let mut bin_buf: Vec<u8> = Vec::new();
        let mut prim_meta: Vec<PrimMeta> = Vec::with_capacity(primitives.len());

        for prim in &primitives {
            let n_verts = prim.positions.len();
            let n_indices = prim.indices.len();

            // Compute bounding box for POSITION accessor.
            let mut pos_min = [f32::INFINITY; 3];
            let mut pos_max = [f32::NEG_INFINITY; 3];
            for p in &prim.positions {
                for i in 0..3 {
                    if p[i] < pos_min[i] {
                        pos_min[i] = p[i];
                    }
                    if p[i] > pos_max[i] {
                        pos_max[i] = p[i];
                    }
                }
            }
            // Fallback for empty primitives.
            if n_verts == 0 {
                pos_min = [0.0; 3];
                pos_max = [0.0; 3];
            }

            let pos_byte_offset = bin_buf.len();
            for p in &prim.positions {
                for &v in p {
                    bin_buf.extend_from_slice(&v.to_le_bytes());
                }
            }

            let norm_byte_offset = bin_buf.len();
            for n in &prim.normals {
                for &v in n {
                    bin_buf.extend_from_slice(&v.to_le_bytes());
                }
            }
            // Pad with zero normals if fewer normals than positions.
            let norm_count = prim.normals.len();
            for _ in norm_count..n_verts {
                bin_buf.extend_from_slice(&[0u8; 12]);
            }

            let uv_byte_offset = bin_buf.len();
            if prim.texcoords.is_empty() {
                // Write zero UVs.
                for _ in 0..n_verts {
                    bin_buf.extend_from_slice(&[0u8; 8]);
                }
            } else {
                for uv in &prim.texcoords {
                    for &v in uv {
                        bin_buf.extend_from_slice(&v.to_le_bytes());
                    }
                }
                // Pad with zeros if fewer UVs than vertices.
                let uv_count = prim.texcoords.len();
                for _ in uv_count..n_verts {
                    bin_buf.extend_from_slice(&[0u8; 8]);
                }
            }

            // Align index buffer to 4-byte boundary.
            let pad_to_align = (4 - bin_buf.len() % 4) % 4;
            bin_buf.extend(std::iter::repeat_n(0u8, pad_to_align));
            let idx_byte_offset = bin_buf.len();
            for &idx in &prim.indices {
                bin_buf.extend_from_slice(&idx.to_le_bytes());
            }

            prim_meta.push(PrimMeta {
                n_verts,
                n_indices,
                pos_byte_offset,
                norm_byte_offset,
                uv_byte_offset,
                idx_byte_offset,
                pos_min,
                pos_max,
            });
        }

        // Pad binary buffer to 4-byte boundary.
        let bin_tail_pad = (4 - bin_buf.len() % 4) % 4;
        bin_buf.extend(std::iter::repeat_n(0u8, bin_tail_pad));
        let total_bin_len = bin_buf.len();

        // Build JSON.
        let json = Self::build_glb_json(scene, &prim_meta, total_bin_len);
        let json_bytes = json.as_bytes();
        let json_raw_len = json_bytes.len();
        let json_padded_len = (json_raw_len + 3) & !3;
        let json_padding = json_padded_len - json_raw_len;

        // GLB: 12-byte header + 8-byte JSON chunk header + JSON + 8-byte BIN chunk header + BIN
        let total_len = 12 + 8 + json_padded_len + 8 + total_bin_len;

        let mut out = Vec::with_capacity(total_len);
        // Header
        out.extend_from_slice(b"glTF");
        out.extend_from_slice(&2u32.to_le_bytes());
        out.extend_from_slice(&(total_len as u32).to_le_bytes());
        // JSON chunk: pad with spaces (0x20) per glTF spec
        out.extend_from_slice(&(json_padded_len as u32).to_le_bytes());
        out.extend_from_slice(&0x4E4F534Au32.to_le_bytes()); // "JSON"
        out.extend_from_slice(json_bytes);
        out.extend(std::iter::repeat_n(0x20u8, json_padding));
        // BIN chunk
        out.extend_from_slice(&(total_bin_len as u32).to_le_bytes());
        out.extend_from_slice(&0x004E4942u32.to_le_bytes()); // "BIN\0"
        out.extend_from_slice(&bin_buf);

        out
    }

    /// Build the glTF JSON for a GLB with a binary buffer.
    fn build_glb_json(scene: &GltfScene, prim_meta: &[PrimMeta], total_bin_len: usize) -> String {
        use std::fmt::Write as FmtWrite;
        let mut out = String::new();

        let _ = writeln!(out, "{{");
        let _ = writeln!(
            out,
            "  \"asset\": {{\"version\": \"2.0\", \"generator\": \"OxiPhysics glTF writer\"}},"
        );
        let _ = writeln!(out, "  \"scene\": 0,");

        // Scenes
        let node_indices: Vec<String> = (0..scene.nodes.len()).map(|i| i.to_string()).collect();
        let _ = writeln!(
            out,
            "  \"scenes\": [{{\"nodes\": [{}]}}],",
            node_indices.join(", ")
        );

        // Nodes
        let _ = writeln!(out, "  \"nodes\": [");
        for (i, node) in scene.nodes.iter().enumerate() {
            let mesh_str = node
                .mesh
                .map(|m| format!(", \"mesh\": {m}"))
                .unwrap_or_default();
            let comma = if i + 1 < scene.nodes.len() { "," } else { "" };
            let _ = writeln!(
                out,
                "    {{\"name\": \"{}\"{mesh_str}}}{comma}",
                escape_json(&node.name)
            );
        }
        let _ = writeln!(out, "  ],");

        // Meshes
        let _ = writeln!(out, "  \"meshes\": [");
        let mut acc_idx = 0usize;
        for (mi, mesh) in scene.meshes.iter().enumerate() {
            let _ = writeln!(
                out,
                "    {{\"name\": \"{}\", \"primitives\": [",
                escape_json(&mesh.name)
            );
            for (pi, _prim) in mesh.primitives.iter().enumerate() {
                let pos_acc = acc_idx;
                let nrm_acc = acc_idx + 1;
                let uv_acc = acc_idx + 2;
                let idx_acc = acc_idx + 3;
                acc_idx += 4;
                let prim_comma = if pi + 1 < mesh.primitives.len() {
                    ","
                } else {
                    ""
                };
                let _ = writeln!(
                    out,
                    "      {{\"attributes\": {{\"POSITION\": {pos_acc}, \"NORMAL\": {nrm_acc}, \"TEXCOORD_0\": {uv_acc}}}, \"indices\": {idx_acc}}}{prim_comma}"
                );
            }
            let mesh_comma = if mi + 1 < scene.meshes.len() { "," } else { "" };
            let _ = writeln!(out, "    ]}}{mesh_comma}");
        }
        let _ = writeln!(out, "  ],");

        // Buffers
        let _ = writeln!(out, "  \"buffers\": [{{\"byteLength\": {total_bin_len}}}],");

        // BufferViews: collect all as strings then join
        let n_meta = prim_meta.len();
        let total_bv = n_meta * 4;
        let mut bv_entries: Vec<String> = Vec::with_capacity(total_bv);
        for meta in prim_meta.iter() {
            let pos_size = meta.n_verts * 12; // [f32;3] = 12 bytes
            let norm_size = meta.n_verts * 12;
            let uv_size = meta.n_verts * 8; // [f32;2] = 8 bytes
            let idx_size = meta.n_indices * 4; // u32 = 4 bytes
            bv_entries.push(format!(
                "    {{\"buffer\": 0, \"byteOffset\": {}, \"byteLength\": {pos_size}, \"target\": 34962}}",
                meta.pos_byte_offset
            ));
            bv_entries.push(format!(
                "    {{\"buffer\": 0, \"byteOffset\": {}, \"byteLength\": {norm_size}, \"target\": 34962}}",
                meta.norm_byte_offset
            ));
            bv_entries.push(format!(
                "    {{\"buffer\": 0, \"byteOffset\": {}, \"byteLength\": {uv_size}, \"target\": 34962}}",
                meta.uv_byte_offset
            ));
            bv_entries.push(format!(
                "    {{\"buffer\": 0, \"byteOffset\": {}, \"byteLength\": {idx_size}, \"target\": 34963}}",
                meta.idx_byte_offset
            ));
        }
        let _ = writeln!(out, "  \"bufferViews\": [");
        let _ = writeln!(out, "{}", bv_entries.join(",\n"));
        let _ = writeln!(out, "  ],");

        // Accessors: collect then join
        let total_acc = n_meta * 4;
        let mut acc_entries: Vec<String> = Vec::with_capacity(total_acc);
        let mut bv_idx = 0usize;
        for meta in prim_meta.iter() {
            let pos_min = meta.pos_min;
            let pos_max = meta.pos_max;
            // POSITION accessor (with min/max bounds)
            acc_entries.push(format!(
                "    {{\"bufferView\": {bv_idx}, \"componentType\": 5126, \"count\": {}, \"type\": \"VEC3\", \
\"min\": [{}, {}, {}], \"max\": [{}, {}, {}]}}",
                meta.n_verts,
                pos_min[0], pos_min[1], pos_min[2],
                pos_max[0], pos_max[1], pos_max[2]
            ));
            bv_idx += 1;
            // NORMAL accessor
            acc_entries.push(format!(
                "    {{\"bufferView\": {bv_idx}, \"componentType\": 5126, \"count\": {}, \"type\": \"VEC3\"}}",
                meta.n_verts
            ));
            bv_idx += 1;
            // TEXCOORD_0 accessor
            acc_entries.push(format!(
                "    {{\"bufferView\": {bv_idx}, \"componentType\": 5126, \"count\": {}, \"type\": \"VEC2\"}}",
                meta.n_verts
            ));
            bv_idx += 1;
            // INDEX accessor
            acc_entries.push(format!(
                "    {{\"bufferView\": {bv_idx}, \"componentType\": 5125, \"count\": {}, \"type\": \"SCALAR\"}}",
                meta.n_indices
            ));
            bv_idx += 1;
        }
        let _ = writeln!(out, "  \"accessors\": [");
        let _ = writeln!(out, "{}", acc_entries.join(",\n"));
        let _ = writeln!(out, "  ]");

        let _ = write!(out, "}}");
        out
    }
}

/// Metadata for a single primitive used during GLB JSON generation.
struct PrimMeta {
    n_verts: usize,
    n_indices: usize,
    pos_byte_offset: usize,
    norm_byte_offset: usize,
    uv_byte_offset: usize,
    idx_byte_offset: usize,
    pos_min: [f32; 3],
    pos_max: [f32; 3],
}
