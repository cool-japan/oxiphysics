//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

use super::functions::compute_normals_from_mesh;
use super::types::{PyColormap, PyLight, PyMaterial, PyParticleRenderer, PySceneNode};
use crate::rasterizer;

/// GPU-ready mesh data with material and optional per-vertex scalar field.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMeshRenderer {
    /// Flat array of vertex positions: \[x0, y0, z0, x1, y1, z1, ...\].
    pub vertices: Vec<f64>,
    /// Flat array of vertex normals (same length as vertices).
    pub normals: Vec<f64>,
    /// Triangle index list (3 indices per triangle).
    pub indices: Vec<u32>,
    /// Surface material.
    pub material: PyMaterial,
    /// Optional per-vertex scalar values for colormap coloring.
    pub scalars: Option<Vec<f64>>,
}
#[pymethods]
impl PyMeshRenderer {
    /// Create a new mesh renderer from raw geometry.
    #[new]
    pub fn new(vertices: Vec<f64>, indices: Vec<u32>, material: PyMaterial) -> Self {
        let normals = compute_normals_from_mesh(vertices.clone(), indices.clone());
        Self {
            vertices,
            normals,
            indices,
            material,
            scalars: None,
        }
    }
    /// Apply a colormap to per-vertex scalar values.
    ///
    /// Sets `material.colormap` and stores the scalar array.
    pub fn apply_colormap(&mut self, scalars: Vec<f64>, colormap: PyColormap) {
        self.scalars = Some(scalars);
        self.material.colormap = Some(colormap);
    }
    /// Returns the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len() / 3
    }
    /// Returns the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
    /// Recompute normals from current vertex/index data.
    pub fn recompute_normals(&mut self) {
        self.normals = compute_normals_from_mesh(self.vertices.clone(), self.indices.clone());
    }
    /// Render to a flat RGBA pixel buffer using a software scanline rasterizer.
    ///
    /// Uses an orthographic projection auto-fitted to the mesh bounding box and
    /// per-vertex Phong shading with a single directional light from `[1,1,1]`.
    /// Returns a buffer of `width * height * 4` bytes (RGBA u8).
    pub fn render_to_buffer(&self, width: u32, height: u32, background: [u8; 4]) -> Vec<u8> {
        let base_color = [
            self.material.base_color[0],
            self.material.base_color[1],
            self.material.base_color[2],
        ];
        let scalars_ref = self.scalars.as_deref();
        if let Some(colormap) = &self.material.colormap {
            let cmap = colormap.clone();
            let map_fn = move |sv: f64| {
                let rgba = cmap.map_scalar(sv);
                [rgba[0], rgba[1], rgba[2]]
            };
            rasterizer::render_mesh(
                rasterizer::MeshRenderData {
                    vertices: &self.vertices,
                    normals: &self.normals,
                    indices: &self.indices,
                    base_color,
                    scalars: scalars_ref,
                    scalar_map_fn: Some(&map_fn),
                },
                width,
                height,
                background,
            )
        } else {
            rasterizer::render_mesh(
                rasterizer::MeshRenderData {
                    vertices: &self.vertices,
                    normals: &self.normals,
                    indices: &self.indices,
                    base_color,
                    scalars: scalars_ref,
                    scalar_map_fn: None,
                },
                width,
                height,
                background,
            )
        }
    }
}
/// Hierarchical scene graph managing nodes, meshes, particles and lights.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySceneGraph {
    /// All nodes in the scene.
    pub nodes: Vec<PySceneNode>,
    /// Mesh renderers attached to node ids.
    pub meshes: Vec<(u32, PyMeshRenderer)>,
    /// Particle renderers attached to node ids.
    pub particles: Vec<(u32, PyParticleRenderer)>,
    /// Lights in the scene.
    pub lights: Vec<PyLight>,
    /// Next free node id.
    next_id: u32,
}
#[pymethods]
impl PySceneGraph {
    /// Create an empty scene graph.
    #[new]
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            meshes: Vec::new(),
            particles: Vec::new(),
            lights: Vec::new(),
            next_id: 0,
        }
    }
    /// Add a new node and return its id.
    pub fn add_node(&mut self, label: String) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.nodes.push(PySceneNode::identity(id, label));
        id
    }
    /// Attach a mesh to a node.
    pub fn add_mesh(&mut self, node_id: u32, mesh: PyMeshRenderer) {
        self.meshes.push((node_id, mesh));
    }
    /// Attach particles to a node.
    pub fn add_particles(&mut self, node_id: u32, particles: PyParticleRenderer) {
        self.particles.push((node_id, particles));
    }
    /// Add a light to the scene.
    pub fn add_light(&mut self, light: PyLight) {
        self.lights.push(light);
    }
    /// Set visibility for a node.
    pub fn set_visible(&mut self, node_id: u32, visible: bool) {
        if let Some(n) = self.nodes.iter_mut().find(|n| n.id == node_id) {
            n.visible = visible;
        }
    }
    /// Set local translation for a node.
    pub fn set_translation(&mut self, node_id: u32, translation: [f64; 3]) {
        if let Some(n) = self.nodes.iter_mut().find(|n| n.id == node_id) {
            n.translation = translation;
        }
    }
    /// Traverse all visible nodes and return their ids in order.
    pub fn traverse_visible(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.visible)
            .map(|n| n.id)
            .collect()
    }
    /// Return total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}
// Default impl is in trait_impls.rs to avoid duplication
