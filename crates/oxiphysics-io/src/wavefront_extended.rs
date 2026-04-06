// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended Wavefront OBJ format with physics metadata.
//!
//! Supports per-vertex velocity annotations and per-material density
//! annotations stored as `# phys_*` comment lines, as well as rigid-body
//! data extraction (centre of mass, inertia tensor) from a mesh with known
//! material densities.

// ── Data types ────────────────────────────────────────────────────────────────

/// A vertex in the extended OBJ mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicsObjVertex {
    /// 3-D position `[x, y, z]`.
    pub pos: [f64; 3],
    /// Optional per-vertex velocity `[vx, vy, vz]`, from `# phys_vel` annotations.
    pub velocity: Option<[f64; 3]>,
}

/// A triangular face referencing three 0-based vertex indices.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicsObjFace {
    /// Vertex indices (0-based).
    pub indices: [usize; 3],
    /// Material name for this face (empty string = no material).
    pub material: String,
}

/// A material entry with an associated density.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicsObjMaterial {
    /// Material name.
    pub name: String,
    /// Mass density in kg/m³.
    pub density: f64,
}

/// Extended OBJ mesh holding geometry and physics metadata.
#[derive(Debug, Clone, Default)]
pub struct PhysicsObjMesh {
    /// Vertex list.
    pub vertices: Vec<PhysicsObjVertex>,
    /// Triangle face list.
    pub faces: Vec<PhysicsObjFace>,
    /// Material table.
    pub materials: Vec<PhysicsObjMaterial>,
}

impl PhysicsObjMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Look up density for a material name, returning `1.0` if not found.
    pub fn density_for(&self, name: &str) -> f64 {
        self.materials
            .iter()
            .find(|m| m.name == name)
            .map(|m| m.density)
            .unwrap_or(1.0)
    }
}

/// Rigid body properties extracted from a mesh.
#[derive(Debug, Clone)]
pub struct RigidBodyData {
    /// Total mass in kg.
    pub mass: f64,
    /// Centre of mass in world space.
    pub com: [f64; 3],
    /// Symmetric inertia tensor (6 components: Ixx, Iyy, Izz, Ixy, Ixz, Iyz).
    pub inertia: [f64; 6],
}

// ── Parsing helpers ───────────────────────────────────────────────────────────

fn parse_f64_triple(parts: &[&str], offset: usize) -> Option<[f64; 3]> {
    if parts.len() < offset + 3 {
        return None;
    }
    let x = parts[offset].parse::<f64>().ok()?;
    let y = parts[offset + 1].parse::<f64>().ok()?;
    let z = parts[offset + 2].parse::<f64>().ok()?;
    Some([x, y, z])
}

fn parse_face_index(s: &str) -> Option<usize> {
    // OBJ face indices can be "v", "v/t", "v/t/n", "v//n" — take the first
    let v_str = s.split('/').next()?;
    let idx: isize = v_str.parse().ok()?;
    if idx > 0 {
        Some((idx - 1) as usize) // 1-based → 0-based
    } else {
        None
    }
}

// ── I/O functions ─────────────────────────────────────────────────────────────

/// Parse an extended Wavefront OBJ string with `# phys_*` annotations.
///
/// Recognised annotations (must appear immediately after the vertex/material
/// they annotate):
/// - `# phys_vel vx vy vz` — velocity for the last declared vertex
/// - `# phys_density name value` — density for material `name`
///
/// Standard OBJ keywords handled: `v`, `f`, `usemtl`.
pub fn read_physics_obj(src: &str) -> PhysicsObjMesh {
    let mut mesh = PhysicsObjMesh::new();
    let mut current_material = String::new();

    for line in src.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();

        match parts[0] {
            "v" => {
                if let Some(pos) = parse_f64_triple(&parts, 1) {
                    mesh.vertices.push(PhysicsObjVertex {
                        pos,
                        velocity: None,
                    });
                }
            }
            "usemtl"
                if parts.len() > 1 => {
                    current_material = parts[1].to_string();
                }
            "f"
                // Only handle triangles (exactly 3 indices)
                if parts.len() >= 4 => {
                    let a = parse_face_index(parts[1]);
                    let b = parse_face_index(parts[2]);
                    let c = parse_face_index(parts[3]);
                    if let (Some(a), Some(b), Some(c)) = (a, b, c) {
                        mesh.faces.push(PhysicsObjFace {
                            indices: [a, b, c],
                            material: current_material.clone(),
                        });
                    }
                }
            "#" if parts.len() >= 2 && parts[1] == "phys_vel" => {
                if let Some(vel) = parse_f64_triple(&parts, 2)
                    && let Some(last) = mesh.vertices.last_mut() {
                        last.velocity = Some(vel);
                    }
            }
            "#" if parts.len() >= 4 && parts[1] == "phys_density" => {
                let mat_name = parts[2].to_string();
                if let Ok(density) = parts[3].parse::<f64>() {
                    // Update existing or push new
                    if let Some(m) = mesh.materials.iter_mut().find(|m| m.name == mat_name) {
                        m.density = density;
                    } else {
                        mesh.materials.push(PhysicsObjMaterial {
                            name: mat_name,
                            density,
                        });
                    }
                }
            }
            _ => {}
        }
    }

    mesh
}

/// Serialise a `PhysicsObjMesh` to the extended OBJ text format.
///
/// Writes:
/// - `v x y z` for each vertex, optionally followed by `# phys_vel vx vy vz`
/// - `# phys_density name value` for each material
/// - `usemtl name` / `f i j k` for each face
pub fn write_physics_obj(mesh: &PhysicsObjMesh) -> String {
    let mut out = String::new();

    // Material density annotations
    for mat in &mesh.materials {
        out.push_str(&material_density_annotation(&mat.name, mat.density));
        out.push('\n');
    }

    // Vertices
    for v in &mesh.vertices {
        out.push_str(&format!("v {} {} {}\n", v.pos[0], v.pos[1], v.pos[2]));
        if let Some(vel) = v.velocity {
            out.push_str(&vertex_velocity_annotation(vel));
            out.push('\n');
        }
    }

    // Faces (grouped by material for readability)
    let mut last_mat = "";
    for face in &mesh.faces {
        if face.material != last_mat {
            if !face.material.is_empty() {
                out.push_str(&format!("usemtl {}\n", face.material));
            }
            last_mat = &face.material;
        }
        // OBJ is 1-based
        out.push_str(&format!(
            "f {} {} {}\n",
            face.indices[0] + 1,
            face.indices[1] + 1,
            face.indices[2] + 1
        ));
    }

    out
}

/// Format a per-vertex velocity as an OBJ comment annotation.
pub fn vertex_velocity_annotation(vel: [f64; 3]) -> String {
    format!("# phys_vel {} {} {}", vel[0], vel[1], vel[2])
}

/// Format a per-material density as an OBJ comment annotation.
pub fn material_density_annotation(name: &str, density: f64) -> String {
    format!("# phys_density {name} {density}")
}

/// Compute rigid-body properties (mass, COM, inertia) from a mesh.
///
/// Each triangle is treated as a thin planar element.  Mass is proportional
/// to triangle area times the density of the face material.  The inertia
/// tensor is approximated from the point masses at triangle centroids.
pub fn extract_rigid_body_data(mesh: &PhysicsObjMesh) -> RigidBodyData {
    let mut total_mass = 0.0_f64;
    let mut com = [0.0_f64; 3];

    // First pass: total mass and COM numerator
    for face in &mesh.faces {
        let [a, b, c] = face.indices;
        if a >= mesh.vertices.len() || b >= mesh.vertices.len() || c >= mesh.vertices.len() {
            continue;
        }
        let pa = mesh.vertices[a].pos;
        let pb = mesh.vertices[b].pos;
        let pc = mesh.vertices[c].pos;

        let area = triangle_area(pa, pb, pc);
        let density = mesh.density_for(&face.material);
        let mass = area * density;
        let centroid = [
            (pa[0] + pb[0] + pc[0]) / 3.0,
            (pa[1] + pb[1] + pc[1]) / 3.0,
            (pa[2] + pb[2] + pc[2]) / 3.0,
        ];

        total_mass += mass;
        com[0] += mass * centroid[0];
        com[1] += mass * centroid[1];
        com[2] += mass * centroid[2];
    }

    if total_mass > 0.0 {
        com[0] /= total_mass;
        com[1] /= total_mass;
        com[2] /= total_mass;
    }

    // Second pass: inertia tensor (Ixx, Iyy, Izz, Ixy, Ixz, Iyz)
    let mut inertia = [0.0_f64; 6];
    for face in &mesh.faces {
        let [a, b, c] = face.indices;
        if a >= mesh.vertices.len() || b >= mesh.vertices.len() || c >= mesh.vertices.len() {
            continue;
        }
        let pa = mesh.vertices[a].pos;
        let pb = mesh.vertices[b].pos;
        let pc = mesh.vertices[c].pos;

        let area = triangle_area(pa, pb, pc);
        let density = mesh.density_for(&face.material);
        let mass = area * density;
        let r = [
            (pa[0] + pb[0] + pc[0]) / 3.0 - com[0],
            (pa[1] + pb[1] + pc[1]) / 3.0 - com[1],
            (pa[2] + pb[2] + pc[2]) / 3.0 - com[2],
        ];

        inertia[0] += mass * (r[1] * r[1] + r[2] * r[2]); // Ixx
        inertia[1] += mass * (r[0] * r[0] + r[2] * r[2]); // Iyy
        inertia[2] += mass * (r[0] * r[0] + r[1] * r[1]); // Izz
        inertia[3] -= mass * r[0] * r[1]; // Ixy
        inertia[4] -= mass * r[0] * r[2]; // Ixz
        inertia[5] -= mass * r[1] * r[2]; // Iyz
    }

    RigidBodyData {
        mass: total_mass,
        com,
        inertia,
    }
}

// ── Internal helper ───────────────────────────────────────────────────────────

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn length3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    0.5 * length3(cross3(ab, ac))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- Annotation helpers ---

    #[test]
    fn test_vertex_velocity_annotation_format() {
        let ann = vertex_velocity_annotation([1.0, 2.0, 3.0]);
        assert!(ann.contains("phys_vel"));
        assert!(ann.contains("1"));
        assert!(ann.contains("2"));
        assert!(ann.contains("3"));
    }

    #[test]
    fn test_material_density_annotation_format() {
        let ann = material_density_annotation("steel", 7800.0);
        assert!(ann.contains("phys_density"));
        assert!(ann.contains("steel"));
        assert!(ann.contains("7800"));
    }

    // --- read_physics_obj ---

    #[test]
    fn test_read_empty_string() {
        let mesh = read_physics_obj("");
        assert!(mesh.vertices.is_empty());
        assert!(mesh.faces.is_empty());
    }

    #[test]
    fn test_read_single_vertex() {
        let src = "v 1.0 2.0 3.0\n";
        let mesh = read_physics_obj(src);
        assert_eq!(mesh.vertices.len(), 1);
        assert_eq!(mesh.vertices[0].pos, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_read_vertex_no_velocity_by_default() {
        let src = "v 0.0 0.0 0.0\n";
        let mesh = read_physics_obj(src);
        assert!(mesh.vertices[0].velocity.is_none());
    }

    #[test]
    fn test_read_phys_vel_annotation() {
        let src = "v 0.0 0.0 0.0\n# phys_vel 1.0 2.0 3.0\n";
        let mesh = read_physics_obj(src);
        assert_eq!(mesh.vertices[0].velocity, Some([1.0, 2.0, 3.0]));
    }

    #[test]
    fn test_read_triangle_face() {
        let src = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let mesh = read_physics_obj(src);
        assert_eq!(mesh.faces.len(), 1);
        assert_eq!(mesh.faces[0].indices, [0, 1, 2]);
    }

    #[test]
    fn test_read_usemtl_sets_material() {
        let src = "v 0 0 0\nv 1 0 0\nv 0 1 0\nusemtl steel\nf 1 2 3\n";
        let mesh = read_physics_obj(src);
        assert_eq!(mesh.faces[0].material, "steel");
    }

    #[test]
    fn test_read_phys_density() {
        let src = "# phys_density wood 500.0\n";
        let mesh = read_physics_obj(src);
        assert_eq!(mesh.materials.len(), 1);
        assert_eq!(mesh.materials[0].name, "wood");
        assert!((mesh.materials[0].density - 500.0).abs() < 1e-9);
    }

    #[test]
    fn test_density_for_unknown_material_returns_one() {
        let mesh = PhysicsObjMesh::new();
        assert!((mesh.density_for("unknown") - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_density_for_known_material() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.materials.push(PhysicsObjMaterial {
            name: "steel".into(),
            density: 7800.0,
        });
        assert!((mesh.density_for("steel") - 7800.0).abs() < 1e-9);
    }

    // --- write_physics_obj ---

    #[test]
    fn test_write_contains_vertex() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.vertices.push(PhysicsObjVertex {
            pos: [1.0, 2.0, 3.0],
            velocity: None,
        });
        let out = write_physics_obj(&mesh);
        assert!(out.contains("v 1"));
    }

    #[test]
    fn test_write_contains_velocity_annotation() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0; 3],
            velocity: Some([4.0, 5.0, 6.0]),
        });
        let out = write_physics_obj(&mesh);
        assert!(out.contains("phys_vel"));
    }

    #[test]
    fn test_write_contains_face() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0; 3],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [1.0, 0.0, 0.0],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 1.0, 0.0],
            velocity: None,
        });
        mesh.faces.push(PhysicsObjFace {
            indices: [0, 1, 2],
            material: String::new(),
        });
        let out = write_physics_obj(&mesh);
        assert!(out.contains("f 1 2 3"));
    }

    #[test]
    fn test_write_contains_density_annotation() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.materials.push(PhysicsObjMaterial {
            name: "iron".into(),
            density: 7874.0,
        });
        let out = write_physics_obj(&mesh);
        assert!(out.contains("phys_density"));
        assert!(out.contains("iron"));
    }

    // --- Roundtrip ---

    fn simple_triangle_mesh() -> PhysicsObjMesh {
        let mut mesh = PhysicsObjMesh::new();
        mesh.materials.push(PhysicsObjMaterial {
            name: "mat".into(),
            density: 1000.0,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 0.0, 0.0],
            velocity: Some([1.0, 0.0, 0.0]),
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [1.0, 0.0, 0.0],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 1.0, 0.0],
            velocity: None,
        });
        mesh.faces.push(PhysicsObjFace {
            indices: [0, 1, 2],
            material: "mat".into(),
        });
        mesh
    }

    #[test]
    fn test_roundtrip_vertex_count() {
        let mesh = simple_triangle_mesh();
        let text = write_physics_obj(&mesh);
        let mesh2 = read_physics_obj(&text);
        assert_eq!(mesh2.vertices.len(), mesh.vertices.len());
    }

    #[test]
    fn test_roundtrip_face_count() {
        let mesh = simple_triangle_mesh();
        let text = write_physics_obj(&mesh);
        let mesh2 = read_physics_obj(&text);
        assert_eq!(mesh2.faces.len(), mesh.faces.len());
    }

    #[test]
    fn test_roundtrip_velocity() {
        let mesh = simple_triangle_mesh();
        let text = write_physics_obj(&mesh);
        let mesh2 = read_physics_obj(&text);
        assert_eq!(mesh2.vertices[0].velocity, Some([1.0, 0.0, 0.0]));
    }

    #[test]
    fn test_roundtrip_material_density() {
        let mesh = simple_triangle_mesh();
        let text = write_physics_obj(&mesh);
        let mesh2 = read_physics_obj(&text);
        assert!((mesh2.density_for("mat") - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_roundtrip_face_indices() {
        let mesh = simple_triangle_mesh();
        let text = write_physics_obj(&mesh);
        let mesh2 = read_physics_obj(&text);
        assert_eq!(mesh2.faces[0].indices, [0, 1, 2]);
    }

    // --- extract_rigid_body_data ---

    #[test]
    fn test_com_single_triangle_at_origin() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 0.0, 0.0],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [2.0, 0.0, 0.0],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 2.0, 0.0],
            velocity: None,
        });
        mesh.faces.push(PhysicsObjFace {
            indices: [0, 1, 2],
            material: String::new(),
        });
        let rb = extract_rigid_body_data(&mesh);
        // Centroid of equilateral-like triangle with these vertices
        assert!(
            (rb.com[0] - 2.0 / 3.0).abs() < 1e-9,
            "com_x = {}",
            rb.com[0]
        );
        assert!(
            (rb.com[1] - 2.0 / 3.0).abs() < 1e-9,
            "com_y = {}",
            rb.com[1]
        );
    }

    #[test]
    fn test_mass_proportional_to_area_and_density() {
        let mut mesh = PhysicsObjMesh::new();
        mesh.materials.push(PhysicsObjMaterial {
            name: "m".into(),
            density: 2.0,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 0.0, 0.0],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [1.0, 0.0, 0.0],
            velocity: None,
        });
        mesh.vertices.push(PhysicsObjVertex {
            pos: [0.0, 1.0, 0.0],
            velocity: None,
        });
        mesh.faces.push(PhysicsObjFace {
            indices: [0, 1, 2],
            material: "m".into(),
        });
        let rb = extract_rigid_body_data(&mesh);
        // Area = 0.5, density = 2.0 => mass = 1.0
        assert!((rb.mass - 1.0).abs() < 1e-9, "mass = {}", rb.mass);
    }

    #[test]
    fn test_empty_mesh_zero_mass() {
        let mesh = PhysicsObjMesh::new();
        let rb = extract_rigid_body_data(&mesh);
        assert!((rb.mass).abs() < 1e-12);
    }

    #[test]
    fn test_inertia_tensor_has_six_components() {
        let mesh = PhysicsObjMesh::new();
        let rb = extract_rigid_body_data(&mesh);
        assert_eq!(rb.inertia.len(), 6);
    }

    #[test]
    fn test_inertia_diagonal_nonnegative() {
        let mesh = simple_triangle_mesh();
        let rb = extract_rigid_body_data(&mesh);
        assert!(rb.inertia[0] >= 0.0, "Ixx negative");
        assert!(rb.inertia[1] >= 0.0, "Iyy negative");
        assert!(rb.inertia[2] >= 0.0, "Izz negative");
    }

    #[test]
    fn test_triangle_area_helper() {
        let a = [0.0; 3];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let area = triangle_area(a, b, c);
        assert!((area - 0.5).abs() < 1e-12, "area = {area}");
    }
}
