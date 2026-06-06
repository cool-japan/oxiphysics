//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Basic IGES file reader.
///
/// Parses the fixed-column IGES format to extract geometric entities.
#[derive(Debug, Clone)]
pub struct IgesReader {
    /// Parsed entities.
    pub entities: Vec<IgesEntity>,
    /// Global section data.
    pub global_params: Vec<String>,
    /// Units identifier.
    pub units: LengthUnit,
}
impl IgesReader {
    /// Create a new IGES reader.
    pub fn new() -> Self {
        Self {
            entities: Vec::new(),
            global_params: Vec::new(),
            units: LengthUnit::Millimeter,
        }
    }
    /// Parse IGES content from a string.
    ///
    /// IGES files have a fixed 80-column format with section identifiers
    /// in column 73.
    pub fn parse(&mut self, content: &str) -> Result<(), String> {
        self.entities.clear();
        self.global_params.clear();
        let lines: Vec<&str> = content.lines().collect();
        let mut _start_lines = Vec::new();
        let mut global_lines = Vec::new();
        let mut directory_lines = Vec::new();
        let mut parameter_lines = Vec::new();
        for line in &lines {
            if line.len() < 73 {
                continue;
            }
            let section = line.chars().nth(72).unwrap_or(' ');
            match section {
                'S' => _start_lines.push(*line),
                'G' => global_lines.push(*line),
                'D' => directory_lines.push(*line),
                'P' => parameter_lines.push(*line),
                _ => {}
            }
        }
        let global_text: String = global_lines
            .iter()
            .map(|l| &l[..72.min(l.len())])
            .collect::<Vec<_>>()
            .join("");
        self.global_params = global_text
            .split(',')
            .map(|s| s.trim().to_string())
            .collect();
        let mut i = 0;
        while i + 1 < directory_lines.len() {
            let line1 = directory_lines[i];
            let _line2 = directory_lines[i + 1];
            let entity_type_str = line1[0..8.min(line1.len())].trim();
            if let Ok(etype_num) = entity_type_str.parse::<u32>() {
                let etype = IgesEntityType::from_type_number(etype_num);
                let seq = i / 2 + 1;
                let entity = IgesEntity::new(seq, etype);
                self.entities.push(entity);
            }
            i += 2;
        }
        let param_text: String = parameter_lines
            .iter()
            .map(|l| &l[..64.min(l.len())])
            .collect::<Vec<_>>()
            .join("");
        let records: Vec<&str> = param_text.split(';').collect();
        for (idx, record) in records.iter().enumerate() {
            if idx >= self.entities.len() {
                break;
            }
            let values: Vec<f64> = record
                .split(',')
                .filter_map(|s| {
                    let s = s.trim();
                    if s.is_empty() {
                        None
                    } else {
                        s.parse::<f64>().ok()
                    }
                })
                .collect();
            self.entities[idx].parameters = values;
        }
        Ok(())
    }
    /// Get all entities of a specific type.
    pub fn find_entities(&self, etype: IgesEntityType) -> Vec<&IgesEntity> {
        self.entities
            .iter()
            .filter(|e| e.entity_type == etype)
            .collect()
    }
    /// Get the number of parsed entities.
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }
    /// Extract points (entity type 116).
    pub fn extract_points(&self) -> Vec<[f64; 3]> {
        let mut points = Vec::new();
        for entity in self.find_entities(IgesEntityType::Point) {
            if entity.parameters.len() >= 3 {
                points.push([
                    entity.parameters[0],
                    entity.parameters[1],
                    entity.parameters[2],
                ]);
            }
        }
        points
    }
    /// Extract lines (entity type 110).
    pub fn extract_lines(&self) -> Vec<([f64; 3], [f64; 3])> {
        let mut lines = Vec::new();
        for entity in self.find_entities(IgesEntityType::Line) {
            if entity.parameters.len() >= 6 {
                let p1 = [
                    entity.parameters[0],
                    entity.parameters[1],
                    entity.parameters[2],
                ];
                let p2 = [
                    entity.parameters[3],
                    entity.parameters[4],
                    entity.parameters[5],
                ];
                lines.push((p1, p2));
            }
        }
        lines
    }
    /// Extract circular arcs (entity type 100).
    pub fn extract_circular_arcs(&self) -> Vec<IgesCircularArc> {
        let mut arcs = Vec::new();
        for entity in self.find_entities(IgesEntityType::CircularArc) {
            if entity.parameters.len() >= 6 {
                arcs.push(IgesCircularArc {
                    z_displacement: entity.parameters[0],
                    center: [entity.parameters[1], entity.parameters[2]],
                    start: [entity.parameters[3], entity.parameters[4]],
                    end: if entity.parameters.len() >= 7 {
                        [entity.parameters[5], entity.parameters[6]]
                    } else {
                        [entity.parameters[3], entity.parameters[4]]
                    },
                });
            }
        }
        arcs
    }
}
/// Circular arc data from IGES entity type 100.
#[derive(Debug, Clone)]
pub struct IgesCircularArc {
    /// Z-displacement of the arc plane.
    pub z_displacement: f64,
    /// Center point \[x, y\].
    pub center: [f64; 2],
    /// Start point \[x, y\].
    pub start: [f64; 2],
    /// End point \[x, y\].
    pub end: [f64; 2],
}
impl IgesCircularArc {
    /// Compute the radius.
    pub fn radius(&self) -> f64 {
        let dx = self.start[0] - self.center[0];
        let dy = self.start[1] - self.center[1];
        (dx * dx + dy * dy).sqrt()
    }
    /// Convert to 3D center point.
    pub fn center_3d(&self) -> [f64; 3] {
        [self.center[0], self.center[1], self.z_displacement]
    }
}
/// A component in an assembly.
#[derive(Debug, Clone)]
pub struct AssemblyComponent {
    /// Component name.
    pub name: String,
    /// BREP solid geometry.
    pub solid: BrepSolid,
    /// Transform relative to parent.
    pub transform: AssemblyTransform,
    /// Child components.
    pub children: Vec<AssemblyComponent>,
}
impl AssemblyComponent {
    /// Create a new component.
    pub fn new(name: &str, solid: BrepSolid, transform: AssemblyTransform) -> Self {
        Self {
            name: name.to_string(),
            solid,
            transform,
            children: Vec::new(),
        }
    }
    /// Add a child component.
    pub fn add_child(&mut self, child: AssemblyComponent) {
        self.children.push(child);
    }
    /// Count total components (including self).
    pub fn total_components(&self) -> usize {
        1 + self
            .children
            .iter()
            .map(|c| c.total_components())
            .sum::<usize>()
    }
    /// Compute the bounding box in world space.
    pub fn world_bounding_box(&self) -> BoundingBox {
        let mut bb = BoundingBox::empty();
        for v in &self.solid.vertices {
            let world_pos = self.transform.apply(v.position);
            bb.include_point(world_pos);
        }
        for child in &self.children {
            let child_bb = child.world_bounding_box();
            bb.include_box(&child_bb);
        }
        bb
    }
    /// Collect all leaf component names.
    pub fn leaf_names(&self) -> Vec<String> {
        if self.children.is_empty() {
            vec![self.name.clone()]
        } else {
            self.children.iter().flat_map(|c| c.leaf_names()).collect()
        }
    }
}
/// Assembly structure holding the root component.
#[derive(Debug, Clone)]
pub struct Assembly {
    /// Assembly name.
    pub name: String,
    /// Root component.
    pub root: AssemblyComponent,
    /// Units used in the assembly.
    pub units: LengthUnit,
}
impl Assembly {
    /// Create a new assembly.
    pub fn new(name: &str, root: AssemblyComponent, units: LengthUnit) -> Self {
        Self {
            name: name.to_string(),
            root,
            units,
        }
    }
    /// Get the total number of components.
    pub fn total_components(&self) -> usize {
        self.root.total_components()
    }
    /// Compute the assembly bounding box.
    pub fn bounding_box(&self) -> BoundingBox {
        self.root.world_bounding_box()
    }
    /// Export the assembly to STL by tessellating all components.
    pub fn to_stl(&self) -> String {
        let exporter = StlExporter::new(&self.name);
        let mut mesh = TriangleMesh::new();
        self.collect_meshes(&self.root, &AssemblyTransform::identity(), &mut mesh);
        exporter.to_ascii_stl(&mesh)
    }
    /// Recursively collect meshes from all components.
    fn collect_meshes(
        &self,
        component: &AssemblyComponent,
        parent_transform: &AssemblyTransform,
        mesh: &mut TriangleMesh,
    ) {
        let world_transform = parent_transform.compose(&component.transform);
        let component_mesh = tessellate_brep(&component.solid, 10);
        let offset = mesh.vertices.len();
        for v in &component_mesh.vertices {
            mesh.vertices.push(world_transform.apply(*v));
        }
        for tri in &component_mesh.triangles {
            mesh.triangles
                .push([tri[0] + offset, tri[1] + offset, tri[2] + offset]);
        }
        for child in &component.children {
            self.collect_meshes(child, &world_transform, mesh);
        }
    }
}
/// Basic STEP file parser for extracting entities.
///
/// Parses ISO 10303-21 (STEP Part 21) format files to extract entity
/// definitions. Does not perform full schema validation.
#[derive(Debug, Clone)]
pub struct StepParser {
    /// Parsed entities.
    pub entities: Vec<StepEntity>,
    /// File description.
    pub description: String,
    /// Schema identifier.
    pub schema: String,
}
impl StepParser {
    /// Create a new empty STEP parser.
    pub fn new() -> Self {
        Self {
            entities: Vec::new(),
            description: String::new(),
            schema: String::new(),
        }
    }
    /// Parse STEP content from a string.
    pub fn parse(&mut self, content: &str) -> Result<(), String> {
        self.entities.clear();
        self.description.clear();
        self.schema.clear();
        let mut in_data = false;
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed == "DATA;" {
                in_data = true;
                continue;
            }
            if trimmed == "ENDSEC;" {
                if in_data {
                    in_data = false;
                }
                continue;
            }
            if trimmed.starts_with("FILE_DESCRIPTION") {
                self.description = trimmed.to_string();
                continue;
            }
            if trimmed.starts_with("FILE_SCHEMA") {
                self.schema = trimmed.to_string();
                continue;
            }
            if in_data
                && trimmed.starts_with('#')
                && let Some(entity) = self.parse_entity_line(trimmed)
            {
                self.entities.push(entity);
            }
        }
        Ok(())
    }
    /// Parse a single entity line like "#123=CARTESIAN_POINT('name',(1.0,2.0,3.0));".
    fn parse_entity_line(&self, line: &str) -> Option<StepEntity> {
        let line = line.trim_end_matches(';');
        let eq_pos = line.find('=')?;
        let id_str = &line[1..eq_pos];
        let id: usize = id_str.parse().ok()?;
        let rest = &line[eq_pos + 1..];
        let paren_pos = rest.find('(');
        let (entity_type, params) = if let Some(pos) = paren_pos {
            let etype = rest[..pos].trim();
            let params = &rest[pos..];
            (etype, params)
        } else {
            (rest.trim(), "")
        };
        Some(StepEntity::new(id, entity_type, params))
    }
    /// Find all entities of a given type.
    pub fn find_entities(&self, entity_type: &str) -> Vec<&StepEntity> {
        self.entities
            .iter()
            .filter(|e| e.entity_type == entity_type)
            .collect()
    }
    /// Find an entity by ID.
    pub fn find_by_id(&self, id: usize) -> Option<&StepEntity> {
        self.entities.iter().find(|e| e.id == id)
    }
    /// Get the total number of parsed entities.
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }
    /// Extract Cartesian points from parsed STEP entities.
    pub fn extract_cartesian_points(&self) -> Vec<(usize, [f64; 3])> {
        let mut points = Vec::new();
        for entity in self.find_entities("CARTESIAN_POINT") {
            if let Some(coords) = self.parse_point_coordinates(&entity.parameters) {
                points.push((entity.id, coords));
            }
        }
        points
    }
    /// Parse point coordinates from a STEP parameter string.
    fn parse_point_coordinates(&self, params: &str) -> Option<[f64; 3]> {
        let inner_start = params.rfind('(')?;
        let inner_end = params[inner_start..].find(')')? + inner_start;
        let coord_str = &params[inner_start + 1..inner_end];
        let coords: Vec<f64> = coord_str
            .split(',')
            .filter_map(|s| s.trim().parse::<f64>().ok())
            .collect();
        if coords.len() >= 3 {
            Some([coords[0], coords[1], coords[2]])
        } else if coords.len() == 2 {
            Some([coords[0], coords[1], 0.0])
        } else {
            None
        }
    }
    /// Extract direction entities.
    pub fn extract_directions(&self) -> Vec<(usize, [f64; 3])> {
        let mut dirs = Vec::new();
        for entity in self.find_entities("DIRECTION") {
            if let Some(dir) = self.parse_point_coordinates(&entity.parameters) {
                dirs.push((entity.id, dir));
            }
        }
        dirs
    }
}
/// Axis-aligned bounding box (AABB) for 3D geometry.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl BoundingBox {
    /// Create an empty (invalid) bounding box.
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }
    /// Create a bounding box from min and max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Expand the bounding box to include a point.
    pub fn include_point(&mut self, p: [f64; 3]) {
        for (i, pi) in p.iter().enumerate() {
            self.min[i] = self.min[i].min(*pi);
            self.max[i] = self.max[i].max(*pi);
        }
    }
    /// Expand to include another bounding box.
    pub fn include_box(&mut self, other: &BoundingBox) {
        for (i, (lo, hi)) in other.min.iter().zip(other.max.iter()).enumerate() {
            self.min[i] = self.min[i].min(*lo);
            self.max[i] = self.max[i].max(*hi);
        }
    }
    /// Compute from a set of points.
    pub fn from_points(points: &[[f64; 3]]) -> Self {
        let mut bb = Self::empty();
        for p in points {
            bb.include_point(*p);
        }
        bb
    }
    /// Get the center of the bounding box.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    /// Get the size (extent) of the bounding box.
    pub fn size(&self) -> [f64; 3] {
        [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ]
    }
    /// Get the diagonal length.
    pub fn diagonal(&self) -> f64 {
        let s = self.size();
        (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt()
    }
    /// Get the volume.
    pub fn volume(&self) -> f64 {
        let s = self.size();
        s[0] * s[1] * s[2]
    }
    /// Get the surface area.
    pub fn surface_area(&self) -> f64 {
        let s = self.size();
        2.0 * (s[0] * s[1] + s[1] * s[2] + s[2] * s[0])
    }
    /// Check if a point is inside the bounding box.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
    /// Check if this box intersects another.
    pub fn intersects(&self, other: &BoundingBox) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Check if the bounding box is valid (non-empty).
    pub fn is_valid(&self) -> bool {
        self.min[0] <= self.max[0] && self.min[1] <= self.max[1] && self.min[2] <= self.max[2]
    }
    /// Apply a unit conversion to the bounding box.
    pub fn convert_units(&self, converter: &UnitConverter) -> BoundingBox {
        BoundingBox {
            min: converter.convert_point(self.min),
            max: converter.convert_point(self.max),
        }
    }
}
/// A vertex in the BREP representation.
#[derive(Debug, Clone)]
pub struct BrepVertex {
    /// Unique identifier.
    pub id: usize,
    /// 3D position.
    pub position: [f64; 3],
    /// Tolerance for vertex coincidence.
    pub tolerance: f64,
}
impl BrepVertex {
    /// Create a new vertex.
    pub fn new(id: usize, position: [f64; 3]) -> Self {
        Self {
            id,
            position,
            tolerance: 1e-6,
        }
    }
}
/// An edge in the BREP representation.
#[derive(Debug, Clone)]
pub struct BrepEdge {
    /// Unique identifier.
    pub id: usize,
    /// Start vertex index.
    pub start_vertex: usize,
    /// End vertex index.
    pub end_vertex: usize,
    /// Curve type.
    pub curve_type: CurveType,
    /// Parameter range \[t_start, t_end\].
    pub parameter_range: [f64; 2],
    /// Intermediate control/sample points.
    pub control_points: Vec<[f64; 3]>,
}
impl BrepEdge {
    /// Create a new linear edge.
    pub fn line(id: usize, start: usize, end: usize) -> Self {
        Self {
            id,
            start_vertex: start,
            end_vertex: end,
            curve_type: CurveType::Line,
            parameter_range: [0.0, 1.0],
            control_points: Vec::new(),
        }
    }
    /// Create a new arc edge.
    pub fn arc(id: usize, start: usize, end: usize, center: [f64; 3]) -> Self {
        Self {
            id,
            start_vertex: start,
            end_vertex: end,
            curve_type: CurveType::CircularArc,
            parameter_range: [0.0, 1.0],
            control_points: vec![center],
        }
    }
    /// Evaluate a point on the edge at parameter t in \[0,1\].
    pub fn evaluate(&self, t: f64, vertices: &[BrepVertex]) -> [f64; 3] {
        let p0 = vertices[self.start_vertex].position;
        let p1 = vertices[self.end_vertex].position;
        match self.curve_type {
            CurveType::Line => lerp3(p0, p1, t),
            CurveType::CircularArc => {
                if let Some(center) = self.control_points.first() {
                    let r0 = sub3(p0, *center);
                    let r1 = sub3(p1, *center);
                    let angle = t * std::f64::consts::PI * 0.5;
                    let cos_a = angle.cos();
                    let sin_a = angle.sin();
                    let r = add3(scale3(r0, cos_a), scale3(r1, sin_a));
                    let radius = (len3(r0) + len3(r1)) * 0.5;
                    let r_len = len3(r);
                    if r_len < 1e-15 {
                        *center
                    } else {
                        add3(*center, scale3(r, radius / r_len))
                    }
                } else {
                    lerp3(p0, p1, t)
                }
            }
            CurveType::BSpline | CurveType::Nurbs => lerp3(p0, p1, t),
        }
    }
    /// Compute the approximate length of the edge.
    pub fn approximate_length(&self, vertices: &[BrepVertex]) -> f64 {
        let n = 20;
        let mut length = 0.0;
        let mut prev = self.evaluate(0.0, vertices);
        for i in 1..=n {
            let t = i as f64 / n as f64;
            let curr = self.evaluate(t, vertices);
            length += len3(sub3(curr, prev));
            prev = curr;
        }
        length
    }
}
/// Surface types in BREP.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceType {
    /// Planar surface.
    Plane,
    /// Cylindrical surface.
    Cylinder,
    /// Spherical surface.
    Sphere,
    /// Conical surface.
    Cone,
    /// Toroidal surface.
    Torus,
    /// B-spline surface.
    BSpline,
}
/// Complete BREP solid.
#[derive(Debug, Clone)]
pub struct BrepSolid {
    /// Name of the solid.
    pub name: String,
    /// Vertices.
    pub vertices: Vec<BrepVertex>,
    /// Edges.
    pub edges: Vec<BrepEdge>,
    /// Faces.
    pub faces: Vec<BrepFace>,
}
impl BrepSolid {
    /// Create a new empty BREP solid.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
        }
    }
    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, position: [f64; 3]) -> usize {
        let id = self.vertices.len();
        self.vertices.push(BrepVertex::new(id, position));
        id
    }
    /// Add a line edge and return its index.
    pub fn add_line_edge(&mut self, start: usize, end: usize) -> usize {
        let id = self.edges.len();
        self.edges.push(BrepEdge::line(id, start, end));
        id
    }
    /// Add a planar face and return its index.
    pub fn add_planar_face(
        &mut self,
        edge_loop: Vec<usize>,
        normal: [f64; 3],
        origin: [f64; 3],
    ) -> usize {
        let id = self.faces.len();
        self.faces
            .push(BrepFace::planar(id, edge_loop, normal, origin));
        id
    }
    /// Compute the bounding box of the solid.
    pub fn bounding_box(&self) -> BoundingBox {
        BoundingBox::from_points(&self.vertices.iter().map(|v| v.position).collect::<Vec<_>>())
    }
    /// Validate the topology: check that all edge vertex references are valid.
    pub fn validate(&self) -> Result<(), String> {
        let n_verts = self.vertices.len();
        let n_edges = self.edges.len();
        for edge in &self.edges {
            if edge.start_vertex >= n_verts {
                return Err(format!(
                    "Edge {} has invalid start vertex {}",
                    edge.id, edge.start_vertex
                ));
            }
            if edge.end_vertex >= n_verts {
                return Err(format!(
                    "Edge {} has invalid end vertex {}",
                    edge.id, edge.end_vertex
                ));
            }
        }
        for face in &self.faces {
            for loop_edges in &face.edge_loops {
                for &eidx in loop_edges {
                    if eidx >= n_edges {
                        return Err(format!("Face {} references invalid edge {}", face.id, eidx));
                    }
                }
            }
        }
        Ok(())
    }
    /// Get Euler characteristic: V - E + F.
    pub fn euler_characteristic(&self) -> i64 {
        self.vertices.len() as i64 - self.edges.len() as i64 + self.faces.len() as i64
    }
    /// Create a box (cuboid) BREP solid.
    pub fn create_box(name: &str, sx: f64, sy: f64, sz: f64) -> Self {
        let mut solid = Self::new(name);
        let v0 = solid.add_vertex([0.0, 0.0, 0.0]);
        let v1 = solid.add_vertex([sx, 0.0, 0.0]);
        let v2 = solid.add_vertex([sx, sy, 0.0]);
        let v3 = solid.add_vertex([0.0, sy, 0.0]);
        let v4 = solid.add_vertex([0.0, 0.0, sz]);
        let v5 = solid.add_vertex([sx, 0.0, sz]);
        let v6 = solid.add_vertex([sx, sy, sz]);
        let v7 = solid.add_vertex([0.0, sy, sz]);
        let e0 = solid.add_line_edge(v0, v1);
        let e1 = solid.add_line_edge(v1, v2);
        let e2 = solid.add_line_edge(v2, v3);
        let e3 = solid.add_line_edge(v3, v0);
        let e4 = solid.add_line_edge(v4, v5);
        let e5 = solid.add_line_edge(v5, v6);
        let e6 = solid.add_line_edge(v6, v7);
        let e7 = solid.add_line_edge(v7, v4);
        let e8 = solid.add_line_edge(v0, v4);
        let e9 = solid.add_line_edge(v1, v5);
        let e10 = solid.add_line_edge(v2, v6);
        let e11 = solid.add_line_edge(v3, v7);
        solid.add_planar_face(vec![e0, e1, e2, e3], [0.0, 0.0, -1.0], [0.0, 0.0, 0.0]);
        solid.add_planar_face(vec![e4, e5, e6, e7], [0.0, 0.0, 1.0], [0.0, 0.0, sz]);
        solid.add_planar_face(vec![e0, e9, e4, e8], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]);
        solid.add_planar_face(vec![e2, e11, e6, e10], [0.0, 1.0, 0.0], [0.0, sy, 0.0]);
        solid.add_planar_face(vec![e3, e8, e7, e11], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        solid.add_planar_face(vec![e1, e10, e5, e9], [1.0, 0.0, 0.0], [sx, 0.0, 0.0]);
        solid
    }
}
/// A face in the BREP representation.
#[derive(Debug, Clone)]
pub struct BrepFace {
    /// Unique identifier.
    pub id: usize,
    /// Edge loops bounding this face (outer loop first, then holes).
    pub edge_loops: Vec<Vec<usize>>,
    /// Surface type.
    pub surface_type: SurfaceType,
    /// Surface normal (for planar faces).
    pub normal: [f64; 3],
    /// Surface origin (for parametric surfaces).
    pub origin: [f64; 3],
    /// Surface parameters (type-dependent).
    pub params: Vec<f64>,
}
impl BrepFace {
    /// Create a new planar face.
    pub fn planar(id: usize, edge_loop: Vec<usize>, normal: [f64; 3], origin: [f64; 3]) -> Self {
        Self {
            id,
            edge_loops: vec![edge_loop],
            surface_type: SurfaceType::Plane,
            normal,
            origin,
            params: Vec::new(),
        }
    }
    /// Create a cylindrical face.
    pub fn cylindrical(
        id: usize,
        edge_loop: Vec<usize>,
        axis: [f64; 3],
        origin: [f64; 3],
        radius: f64,
    ) -> Self {
        Self {
            id,
            edge_loops: vec![edge_loop],
            surface_type: SurfaceType::Cylinder,
            normal: axis,
            origin,
            params: vec![radius],
        }
    }
    /// Create a spherical face.
    pub fn spherical(id: usize, edge_loop: Vec<usize>, center: [f64; 3], radius: f64) -> Self {
        Self {
            id,
            edge_loops: vec![edge_loop],
            surface_type: SurfaceType::Sphere,
            normal: [0.0, 0.0, 1.0],
            origin: center,
            params: vec![radius],
        }
    }
    /// Get the outer edge loop.
    pub fn outer_loop(&self) -> &[usize] {
        &self.edge_loops[0]
    }
    /// Get hole loops (if any).
    pub fn hole_loops(&self) -> &[Vec<usize>] {
        if self.edge_loops.len() > 1 {
            &self.edge_loops[1..]
        } else {
            &[]
        }
    }
}
/// IGES entity types (100-199 range: geometry).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IgesEntityType {
    /// Type 100: Circular Arc.
    CircularArc,
    /// Type 102: Composite Curve.
    CompositeCurve,
    /// Type 104: Conic Arc.
    ConicArc,
    /// Type 106: Copious Data.
    CopiousData,
    /// Type 108: Plane.
    Plane,
    /// Type 110: Line.
    Line,
    /// Type 112: Parametric Spline Curve.
    ParametricSpline,
    /// Type 114: Parametric Spline Surface.
    ParametricSplineSurface,
    /// Type 116: Point.
    Point,
    /// Type 118: Ruled Surface.
    RuledSurface,
    /// Type 120: Surface of Revolution.
    SurfaceOfRevolution,
    /// Type 122: Tabulated Cylinder.
    TabulatedCylinder,
    /// Type 124: Transformation Matrix.
    TransformationMatrix,
    /// Type 126: Rational B-Spline Curve.
    RationalBSplineCurve,
    /// Type 128: Rational B-Spline Surface.
    RationalBSplineSurface,
    /// Type 130: Offset Curve.
    OffsetCurve,
    /// Type 140: Offset Surface.
    OffsetSurface,
    /// Type 142: Curve on Parametric Surface.
    CurveOnSurface,
    /// Type 144: Trimmed Parametric Surface.
    TrimmedSurface,
    /// Unknown entity type.
    Unknown(u32),
}
impl IgesEntityType {
    /// Convert from entity type number.
    pub fn from_type_number(num: u32) -> Self {
        match num {
            100 => IgesEntityType::CircularArc,
            102 => IgesEntityType::CompositeCurve,
            104 => IgesEntityType::ConicArc,
            106 => IgesEntityType::CopiousData,
            108 => IgesEntityType::Plane,
            110 => IgesEntityType::Line,
            112 => IgesEntityType::ParametricSpline,
            114 => IgesEntityType::ParametricSplineSurface,
            116 => IgesEntityType::Point,
            118 => IgesEntityType::RuledSurface,
            120 => IgesEntityType::SurfaceOfRevolution,
            122 => IgesEntityType::TabulatedCylinder,
            124 => IgesEntityType::TransformationMatrix,
            126 => IgesEntityType::RationalBSplineCurve,
            128 => IgesEntityType::RationalBSplineSurface,
            130 => IgesEntityType::OffsetCurve,
            140 => IgesEntityType::OffsetSurface,
            142 => IgesEntityType::CurveOnSurface,
            144 => IgesEntityType::TrimmedSurface,
            _ => IgesEntityType::Unknown(num),
        }
    }
    /// Convert to entity type number.
    pub fn to_type_number(&self) -> u32 {
        match self {
            IgesEntityType::CircularArc => 100,
            IgesEntityType::CompositeCurve => 102,
            IgesEntityType::ConicArc => 104,
            IgesEntityType::CopiousData => 106,
            IgesEntityType::Plane => 108,
            IgesEntityType::Line => 110,
            IgesEntityType::ParametricSpline => 112,
            IgesEntityType::ParametricSplineSurface => 114,
            IgesEntityType::Point => 116,
            IgesEntityType::RuledSurface => 118,
            IgesEntityType::SurfaceOfRevolution => 120,
            IgesEntityType::TabulatedCylinder => 122,
            IgesEntityType::TransformationMatrix => 124,
            IgesEntityType::RationalBSplineCurve => 126,
            IgesEntityType::RationalBSplineSurface => 128,
            IgesEntityType::OffsetCurve => 130,
            IgesEntityType::OffsetSurface => 140,
            IgesEntityType::CurveOnSurface => 142,
            IgesEntityType::TrimmedSurface => 144,
            IgesEntityType::Unknown(n) => *n,
        }
    }
    /// Check if this is a curve entity.
    pub fn is_curve(&self) -> bool {
        matches!(
            self,
            IgesEntityType::CircularArc
                | IgesEntityType::CompositeCurve
                | IgesEntityType::ConicArc
                | IgesEntityType::Line
                | IgesEntityType::ParametricSpline
                | IgesEntityType::RationalBSplineCurve
                | IgesEntityType::OffsetCurve
        )
    }
    /// Check if this is a surface entity.
    pub fn is_surface(&self) -> bool {
        matches!(
            self,
            IgesEntityType::Plane
                | IgesEntityType::RuledSurface
                | IgesEntityType::SurfaceOfRevolution
                | IgesEntityType::TabulatedCylinder
                | IgesEntityType::RationalBSplineSurface
                | IgesEntityType::ParametricSplineSurface
                | IgesEntityType::OffsetSurface
                | IgesEntityType::TrimmedSurface
        )
    }
}
/// Length units supported by CAD files.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LengthUnit {
    /// Meters (SI base unit).
    Meter,
    /// Millimeters.
    Millimeter,
    /// Centimeters.
    Centimeter,
    /// Inches.
    Inch,
    /// Feet.
    Foot,
    /// Micrometers.
    Micrometer,
}
impl LengthUnit {
    /// Conversion factor to meters.
    pub fn to_meters_factor(&self) -> f64 {
        match self {
            LengthUnit::Meter => 1.0,
            LengthUnit::Millimeter => 0.001,
            LengthUnit::Centimeter => 0.01,
            LengthUnit::Inch => 0.0254,
            LengthUnit::Foot => 0.3048,
            LengthUnit::Micrometer => 1e-6,
        }
    }
    /// Convert a length from this unit to another unit.
    pub fn convert(&self, value: f64, to: LengthUnit) -> f64 {
        value * self.to_meters_factor() / to.to_meters_factor()
    }
    /// Convert a 3D point from this unit to another.
    pub fn convert_point(&self, point: [f64; 3], to: LengthUnit) -> [f64; 3] {
        let factor = self.to_meters_factor() / to.to_meters_factor();
        [point[0] * factor, point[1] * factor, point[2] * factor]
    }
}
/// Unit converter for CAD data.
#[derive(Debug, Clone)]
pub struct UnitConverter {
    /// Source unit.
    pub from: LengthUnit,
    /// Target unit.
    pub to: LengthUnit,
}
impl UnitConverter {
    /// Create a new unit converter.
    pub fn new(from: LengthUnit, to: LengthUnit) -> Self {
        Self { from, to }
    }
    /// Get the conversion factor.
    pub fn factor(&self) -> f64 {
        self.from.to_meters_factor() / self.to.to_meters_factor()
    }
    /// Convert a scalar value.
    pub fn convert_scalar(&self, value: f64) -> f64 {
        value * self.factor()
    }
    /// Convert a 3D point.
    pub fn convert_point(&self, point: [f64; 3]) -> [f64; 3] {
        let f = self.factor();
        [point[0] * f, point[1] * f, point[2] * f]
    }
    /// Convert an array of points.
    pub fn convert_points(&self, points: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let f = self.factor();
        points
            .iter()
            .map(|p| [p[0] * f, p[1] * f, p[2] * f])
            .collect()
    }
}
/// STL file writer for triangle meshes.
#[derive(Debug, Clone)]
pub struct StlExporter {
    /// Name of the solid in the STL file.
    pub solid_name: String,
}
impl StlExporter {
    /// Create a new STL exporter.
    pub fn new(solid_name: &str) -> Self {
        Self {
            solid_name: solid_name.to_string(),
        }
    }
    /// Generate ASCII STL content from a triangle mesh.
    pub fn to_ascii_stl(&self, mesh: &TriangleMesh) -> String {
        let mut output = format!("solid {}\n", self.solid_name);
        for (i, tri) in mesh.triangles.iter().enumerate() {
            let p0 = mesh.vertices[tri[0]];
            let p1 = mesh.vertices[tri[1]];
            let p2 = mesh.vertices[tri[2]];
            let normal = if i < mesh.normals.len() {
                mesh.normals[i]
            } else {
                let e1 = sub3(p1, p0);
                let e2 = sub3(p2, p0);
                normalize3(cross3(e1, e2))
            };
            output += &format!(
                "  facet normal {:.6} {:.6} {:.6}\n",
                normal[0], normal[1], normal[2]
            );
            output += "    outer loop\n";
            output += &format!("      vertex {:.6} {:.6} {:.6}\n", p0[0], p0[1], p0[2]);
            output += &format!("      vertex {:.6} {:.6} {:.6}\n", p1[0], p1[1], p1[2]);
            output += &format!("      vertex {:.6} {:.6} {:.6}\n", p2[0], p2[1], p2[2]);
            output += "    endloop\n";
            output += "  endfacet\n";
        }
        output += &format!("endsolid {}\n", self.solid_name);
        output
    }
    /// Generate binary STL content as bytes.
    pub fn to_binary_stl(&self, mesh: &TriangleMesh) -> Vec<u8> {
        let mut data = Vec::new();
        let mut header = [0u8; 80];
        let name_bytes = self.solid_name.as_bytes();
        let copy_len = name_bytes.len().min(80);
        header[..copy_len].copy_from_slice(&name_bytes[..copy_len]);
        data.extend_from_slice(&header);
        let n_tri = mesh.triangles.len() as u32;
        data.extend_from_slice(&n_tri.to_le_bytes());
        for (i, tri) in mesh.triangles.iter().enumerate() {
            let p0 = mesh.vertices[tri[0]];
            let p1 = mesh.vertices[tri[1]];
            let p2 = mesh.vertices[tri[2]];
            let normal = if i < mesh.normals.len() {
                mesh.normals[i]
            } else {
                let e1 = sub3(p1, p0);
                let e2 = sub3(p2, p0);
                normalize3(cross3(e1, e2))
            };
            for &c in &normal {
                data.extend_from_slice(&(c as f32).to_le_bytes());
            }
            for p in &[p0, p1, p2] {
                for &c in p {
                    data.extend_from_slice(&(c as f32).to_le_bytes());
                }
            }
            data.extend_from_slice(&0u16.to_le_bytes());
        }
        data
    }
    /// Export a BREP solid as ASCII STL.
    pub fn export_brep_ascii(&self, solid: &BrepSolid) -> String {
        let mesh = tessellate_brep(solid, 10);
        self.to_ascii_stl(&mesh)
    }
}
/// A parsed STEP entity.
#[derive(Debug, Clone)]
pub struct StepEntity {
    /// Entity ID (e.g., #123).
    pub id: usize,
    /// Entity type name (e.g., "CARTESIAN_POINT").
    pub entity_type: String,
    /// Raw parameter string.
    pub parameters: String,
}
impl StepEntity {
    /// Create a new STEP entity.
    pub fn new(id: usize, entity_type: &str, parameters: &str) -> Self {
        Self {
            id,
            entity_type: entity_type.to_string(),
            parameters: parameters.to_string(),
        }
    }
}
/// Parsed IGES entity.
#[derive(Debug, Clone)]
pub struct IgesEntity {
    /// Directory entry sequence number.
    pub sequence_number: usize,
    /// Entity type.
    pub entity_type: IgesEntityType,
    /// Parameter data pointer.
    pub parameter_pointer: usize,
    /// Status flags.
    pub status: u32,
    /// Entity label.
    pub label: String,
    /// Raw parameter data values.
    pub parameters: Vec<f64>,
}
impl IgesEntity {
    /// Create a new IGES entity.
    pub fn new(seq: usize, etype: IgesEntityType) -> Self {
        Self {
            sequence_number: seq,
            entity_type: etype,
            parameter_pointer: 0,
            status: 0,
            label: String::new(),
            parameters: Vec::new(),
        }
    }
}
/// Transform for positioning a component in an assembly.
#[derive(Debug, Clone)]
pub struct AssemblyTransform {
    /// 3x3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Translation vector.
    pub translation: [f64; 3],
}
impl AssemblyTransform {
    /// Create an identity transform.
    pub fn identity() -> Self {
        Self {
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
        }
    }
    /// Create a translation-only transform.
    pub fn translation(t: [f64; 3]) -> Self {
        Self {
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: t,
        }
    }
    /// Apply the transform to a point.
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let rotated = [
            self.rotation[0][0] * p[0] + self.rotation[0][1] * p[1] + self.rotation[0][2] * p[2],
            self.rotation[1][0] * p[0] + self.rotation[1][1] * p[1] + self.rotation[1][2] * p[2],
            self.rotation[2][0] * p[0] + self.rotation[2][1] * p[1] + self.rotation[2][2] * p[2],
        ];
        add3(rotated, self.translation)
    }
    /// Compose two transforms: self * other.
    pub fn compose(&self, other: &AssemblyTransform) -> AssemblyTransform {
        let mut new_rot = [[0.0; 3]; 3];
        for (i, row) in new_rot.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                for (k, &self_ik) in self.rotation[i].iter().enumerate() {
                    *cell += self_ik * other.rotation[k][j];
                }
            }
        }
        let new_trans = self.apply(other.translation);
        AssemblyTransform {
            rotation: new_rot,
            translation: new_trans,
        }
    }
    /// Get the inverse transform.
    pub fn inverse(&self) -> AssemblyTransform {
        let r_inv = [
            [
                self.rotation[0][0],
                self.rotation[1][0],
                self.rotation[2][0],
            ],
            [
                self.rotation[0][1],
                self.rotation[1][1],
                self.rotation[2][1],
            ],
            [
                self.rotation[0][2],
                self.rotation[1][2],
                self.rotation[2][2],
            ],
        ];
        let t_inv = [
            -(r_inv[0][0] * self.translation[0]
                + r_inv[0][1] * self.translation[1]
                + r_inv[0][2] * self.translation[2]),
            -(r_inv[1][0] * self.translation[0]
                + r_inv[1][1] * self.translation[1]
                + r_inv[1][2] * self.translation[2]),
            -(r_inv[2][0] * self.translation[0]
                + r_inv[2][1] * self.translation[1]
                + r_inv[2][2] * self.translation[2]),
        ];
        AssemblyTransform {
            rotation: r_inv,
            translation: t_inv,
        }
    }
}
/// Triangle mesh from tessellation.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices (each triple is one triangle).
    pub triangles: Vec<[usize; 3]>,
    /// Per-triangle normals (optional).
    pub normals: Vec<[f64; 3]>,
}
impl TriangleMesh {
    /// Create a new empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
            normals: Vec::new(),
        }
    }
    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, position: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(position);
        idx
    }
    /// Add a triangle.
    pub fn add_triangle(&mut self, v0: usize, v1: usize, v2: usize) {
        self.triangles.push([v0, v1, v2]);
    }
    /// Compute normals for all triangles.
    pub fn compute_normals(&mut self) {
        self.normals.clear();
        for tri in &self.triangles {
            let p0 = self.vertices[tri[0]];
            let p1 = self.vertices[tri[1]];
            let p2 = self.vertices[tri[2]];
            let e1 = sub3(p1, p0);
            let e2 = sub3(p2, p0);
            let n = normalize3(cross3(e1, e2));
            self.normals.push(n);
        }
    }
    /// Compute the bounding box.
    pub fn bounding_box(&self) -> BoundingBox {
        BoundingBox::from_points(&self.vertices)
    }
    /// Compute the total surface area.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0;
        for tri in &self.triangles {
            let p0 = self.vertices[tri[0]];
            let p1 = self.vertices[tri[1]];
            let p2 = self.vertices[tri[2]];
            let e1 = sub3(p1, p0);
            let e2 = sub3(p2, p0);
            area += len3(cross3(e1, e2)) * 0.5;
        }
        area
    }
    /// Count the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    /// Count the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }
    /// Merge another mesh into this one.
    pub fn merge(&mut self, other: &TriangleMesh) {
        let offset = self.vertices.len();
        self.vertices.extend_from_slice(&other.vertices);
        for tri in &other.triangles {
            self.triangles
                .push([tri[0] + offset, tri[1] + offset, tri[2] + offset]);
        }
        self.normals.extend_from_slice(&other.normals);
    }
}
/// Types of curves in BREP.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurveType {
    /// Straight line segment.
    Line,
    /// Circular arc.
    CircularArc,
    /// B-spline curve.
    BSpline,
    /// NURBS curve.
    Nurbs,
}
