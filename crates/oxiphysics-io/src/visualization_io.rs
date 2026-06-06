// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Visualization output formats for physics simulations.
//!
//! Supports: ParaView state files (.pvsm), VisIt database, Blender Python export,
//! Matplotlib JSON descriptor, D3.js data export, WebGL buffer export,
//! glTF 2.0 physics annotations, VDB sparse volume, OpenEXR multi-channel,
//! and Cinema database format.

use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

// ─────────────────────────────────────────────────────────────────────────────
// Shared geometry types (no nalgebra — use [f64;3] arrays)
// ─────────────────────────────────────────────────────────────────────────────

/// A 3D point using plain array storage.
#[derive(Clone, Debug, PartialEq)]
pub struct Point3 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
}

impl Point3 {
    /// Create a new [`Point3`].
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Convert to array.
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Distance to another point.
    pub fn dist(&self, other: &Point3) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// A scalar field on a regular Cartesian grid.
#[derive(Clone, Debug)]
pub struct ScalarField3D {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Number of cells in z.
    pub nz: usize,
    /// Cell spacing.
    pub dx: f64,
    /// Origin.
    pub origin: [f64; 3],
    /// Flat (x-major) data array.
    pub data: Vec<f64>,
}

impl ScalarField3D {
    /// Create a zero-initialised [`ScalarField3D`].
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, origin: [f64; 3]) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            data: vec![0.0; nx * ny * nz],
        }
    }

    /// Index into flat array.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.ny * self.nz + iy * self.nz + iz
    }

    /// Set a value.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, val: f64) {
        let i = self.idx(ix, iy, iz);
        self.data[i] = val;
    }

    /// Get a value.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.data[self.idx(ix, iy, iz)]
    }

    /// Minimum value.
    pub fn min_val(&self) -> f64 {
        self.data.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Maximum value.
    pub fn max_val(&self) -> f64 {
        self.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ParaView State File (.pvsm) writer
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for a ParaView state file.
#[derive(Clone, Debug)]
pub struct ParaviewStateConfig {
    /// Source VTU file path.
    pub source_file: String,
    /// Name of the scalar array to colour by.
    pub colour_array: String,
    /// Colour map preset name (e.g., "Rainbow Desaturated").
    pub colormap: String,
    /// Scalar range \[min, max\].
    pub scalar_range: [f64; 2],
    /// Camera position.
    pub camera_pos: [f64; 3],
    /// Camera focal point.
    pub camera_focal: [f64; 3],
    /// Camera view-up vector.
    pub camera_up: [f64; 3],
    /// Image width in pixels.
    pub image_width: u32,
    /// Image height in pixels.
    pub image_height: u32,
}

impl ParaviewStateConfig {
    /// Create a default [`ParaviewStateConfig`].
    pub fn default_config(source_file: impl Into<String>) -> Self {
        Self {
            source_file: source_file.into(),
            colour_array: "pressure".into(),
            colormap: "Cool to Warm".into(),
            scalar_range: [0.0, 1.0],
            camera_pos: [0.0, 0.0, 5.0],
            camera_focal: [0.0, 0.0, 0.0],
            camera_up: [0.0, 1.0, 0.0],
            image_width: 1920,
            image_height: 1080,
        }
    }
}

/// Write a minimal ParaView state (.pvsm) XML file to a string.
pub fn write_paraview_state(cfg: &ParaviewStateConfig) -> String {
    let mut out = String::new();
    let _ = writeln!(out, r#"<?xml version="1.0"?>"#);
    let _ = writeln!(out, r#"<ParaView>"#);
    let _ = writeln!(out, r#"  <ServerManagerState version="5.10.0">"#);
    // Source
    let _ = writeln!(
        out,
        r#"    <Proxy group="sources" type="XMLUnstructuredGridReader" id="1001">"#
    );
    let _ = writeln!(
        out,
        r#"      <Property name="FileName" number_of_elements="1">"#
    );
    let _ = writeln!(
        out,
        r#"        <Element index="0" value="{}"/>"#,
        cfg.source_file
    );
    let _ = writeln!(out, r#"      </Property>"#);
    let _ = writeln!(out, r#"    </Proxy>"#);
    // Display
    let _ = writeln!(
        out,
        r#"    <Proxy group="representations" type="GeometryRepresentation" id="2001">"#
    );
    let _ = writeln!(out, r#"      <Property name="ColorArrayName">"#);
    let _ = writeln!(out, r#"        <Element index="0" value="POINTS"/>"#);
    let _ = writeln!(
        out,
        r#"        <Element index="1" value="{}"/>"#,
        cfg.colour_array
    );
    let _ = writeln!(out, r#"      </Property>"#);
    let _ = writeln!(out, r#"      <Property name="LookupTable" value="3001"/>"#);
    let _ = writeln!(out, r#"    </Proxy>"#);
    // LUT
    let _ = writeln!(
        out,
        r#"    <Proxy group="lookup_tables" type="PVLookupTable" id="3001">"#
    );
    let _ = writeln!(
        out,
        r#"      <Property name="ColorSpace"><Element index="0" value="HSV"/></Property>"#
    );
    let _ = writeln!(
        out,
        r#"      <Property name="ScalarRangeInitialized"><Element index="0" value="1"/></Property>"#
    );
    let _ = writeln!(
        out,
        r#"      <Property name="RGBPoints" number_of_elements="8">"#
    );
    let _ = writeln!(
        out,
        r#"        <Element index="0" value="{}"/><Element index="1" value="0.23137"/>"#,
        cfg.scalar_range[0]
    );
    let _ = writeln!(
        out,
        r#"        <Element index="2" value="0.29803"/><Element index="3" value="0.75294"/>"#
    );
    let _ = writeln!(
        out,
        r#"        <Element index="4" value="{}"/><Element index="5" value="0.70588"/>"#,
        cfg.scalar_range[1]
    );
    let _ = writeln!(
        out,
        r#"        <Element index="6" value="0.01568"/><Element index="7" value="0.14901"/>"#
    );
    let _ = writeln!(out, r#"      </Property>"#);
    let _ = writeln!(out, r#"    </Proxy>"#);
    // Camera
    let _ = writeln!(
        out,
        r#"    <Proxy group="views" type="RenderView" id="4001">"#
    );
    let _ = writeln!(
        out,
        r#"      <Property name="CameraPosition" number_of_elements="3"><Element index="0" value="{}"/><Element index="1" value="{}"/><Element index="2" value="{}"/></Property>"#,
        cfg.camera_pos[0], cfg.camera_pos[1], cfg.camera_pos[2]
    );
    let _ = writeln!(
        out,
        r#"      <Property name="CameraFocalPoint" number_of_elements="3"><Element index="0" value="{}"/><Element index="1" value="{}"/><Element index="2" value="{}"/></Property>"#,
        cfg.camera_focal[0], cfg.camera_focal[1], cfg.camera_focal[2]
    );
    let _ = writeln!(
        out,
        r#"      <Property name="ViewSize" number_of_elements="2"><Element index="0" value="{}"/><Element index="1" value="{}"/></Property>"#,
        cfg.image_width, cfg.image_height
    );
    let _ = writeln!(out, r#"    </Proxy>"#);
    let _ = writeln!(out, r#"  </ServerManagerState>"#);
    let _ = writeln!(out, r#"</ParaView>"#);
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// VisIt database format writer
// ─────────────────────────────────────────────────────────────────────────────

/// A VisIt database descriptor for a time series.
#[derive(Clone, Debug)]
pub struct VisItDatabase {
    /// Base name of the simulation files.
    pub basename: String,
    /// File extension (e.g., "vtk").
    pub extension: String,
    /// Time steps.
    pub times: Vec<f64>,
    /// Variable names present in each file.
    pub variables: Vec<String>,
    /// Spatial dimension (2 or 3).
    pub ndim: u8,
}

impl VisItDatabase {
    /// Create a new [`VisItDatabase`] descriptor.
    pub fn new(basename: impl Into<String>, extension: impl Into<String>, ndim: u8) -> Self {
        Self {
            basename: basename.into(),
            extension: extension.into(),
            times: Vec::new(),
            variables: Vec::new(),
            ndim,
        }
    }

    /// Add a time step.
    pub fn add_time(&mut self, t: f64) {
        self.times.push(t);
    }

    /// Add a variable name.
    pub fn add_variable(&mut self, name: impl Into<String>) {
        self.variables.push(name.into());
    }

    /// Write a .visit index file for VisIt.
    pub fn write_visit_index(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "!NBLOCKS 1");
        for (i, _t) in self.times.iter().enumerate() {
            let _ = writeln!(out, "{}{:06}.{}", self.basename, i, self.extension);
        }
        out
    }

    /// Write a summary manifest (JSON-like) for the database.
    pub fn write_manifest(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "{{");
        let _ = writeln!(out, r#"  "basename": "{}","#, self.basename);
        let _ = writeln!(out, r#"  "ndim": {},"#, self.ndim);
        let _ = writeln!(out, r#"  "n_times": {},"#, self.times.len());
        let _ = writeln!(out, r#"  "variables": ["#);
        for (i, v) in self.variables.iter().enumerate() {
            let comma = if i + 1 < self.variables.len() {
                ","
            } else {
                ""
            };
            let _ = writeln!(out, r#"    "{}"{}"#, v, comma);
        }
        let _ = writeln!(out, r#"  ],"#);
        let _ = writeln!(out, r#"  "times": ["#);
        for (i, t) in self.times.iter().enumerate() {
            let comma = if i + 1 < self.times.len() { "," } else { "" };
            let _ = writeln!(out, "    {}{}", t, comma);
        }
        let _ = writeln!(out, r#"  ]"#);
        let _ = writeln!(out, "}}");
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Blender Python export
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for a Blender bpy export script.
#[derive(Clone, Debug)]
pub struct BlenderExportConfig {
    /// Name of the Blender object.
    pub object_name: String,
    /// Output .blend file path.
    pub output_path: String,
    /// Whether to add a rigid-body physics modifier.
    pub rigid_body: bool,
    /// Whether to add subdivision surface.
    pub subdivision: bool,
    /// Subdivision level.
    pub subdiv_level: u32,
    /// Background colour (RGBA).
    pub bg_colour: [f64; 4],
    /// Emission strength for volume shaders.
    pub emission_strength: f64,
}

impl BlenderExportConfig {
    /// Create a default [`BlenderExportConfig`].
    pub fn new(object_name: impl Into<String>, output_path: impl Into<String>) -> Self {
        Self {
            object_name: object_name.into(),
            output_path: output_path.into(),
            rigid_body: false,
            subdivision: false,
            subdiv_level: 2,
            bg_colour: [0.05, 0.05, 0.05, 1.0],
            emission_strength: 1.0,
        }
    }

    /// Generate a Blender bpy Python script string.
    pub fn generate_script(&self, mesh_file: &str) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "import bpy");
        let _ = writeln!(s, "import bmesh");
        let _ = writeln!(s);
        let _ = writeln!(s, "# Clear existing objects");
        let _ = writeln!(s, "bpy.ops.object.select_all(action='SELECT')");
        let _ = writeln!(s, "bpy.ops.object.delete()");
        let _ = writeln!(s);
        let _ = writeln!(s, "# Import mesh");
        let _ = writeln!(s, "bpy.ops.import_scene.obj(filepath='{}')", mesh_file);
        let _ = writeln!(s, "obj = bpy.context.selected_objects[0]");
        let _ = writeln!(s, "obj.name = '{}'", self.object_name);
        if self.rigid_body {
            let _ = writeln!(s, "bpy.ops.rigidbody.object_add()");
            let _ = writeln!(s, "obj.rigid_body.type = 'ACTIVE'");
        }
        if self.subdivision {
            let _ = writeln!(s, "mod = obj.modifiers.new(name='Subdiv', type='SUBSURF')");
            let _ = writeln!(s, "mod.levels = {}", self.subdiv_level);
            let _ = writeln!(s, "mod.render_levels = {}", self.subdiv_level);
        }
        let _ = writeln!(s, "# World background");
        let _ = writeln!(
            s,
            "bpy.context.scene.world.node_tree.nodes['Background'].inputs[0].default_value = ({}, {}, {}, {})",
            self.bg_colour[0], self.bg_colour[1], self.bg_colour[2], self.bg_colour[3]
        );
        let _ = writeln!(s, "# Save blend file");
        let _ = writeln!(
            s,
            "bpy.ops.wm.save_as_mainfile(filepath='{}')",
            self.output_path
        );
        s
    }

    /// Generate a Python script for volume rendering.
    pub fn generate_volume_script(&self, vdb_file: &str) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "import bpy");
        let _ = writeln!(s, "bpy.ops.object.volume_import(filepath='{}')", vdb_file);
        let _ = writeln!(s, "vol = bpy.context.active_object");
        let _ = writeln!(s, "mat = bpy.data.materials.new(name='VolMat')");
        let _ = writeln!(s, "mat.use_nodes = True");
        let _ = writeln!(s, "nodes = mat.node_tree.nodes");
        let _ = writeln!(s, "nodes.clear()");
        let _ = writeln!(s, "emission = nodes.new('ShaderNodeEmission')");
        let _ = writeln!(
            s,
            "emission.inputs['Strength'].default_value = {}",
            self.emission_strength
        );
        let _ = writeln!(s, "vol.data.materials.append(mat)");
        s
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Matplotlib JSON descriptor
// ─────────────────────────────────────────────────────────────────────────────

/// A descriptor for a Matplotlib figure, exportable as JSON.
#[derive(Clone, Debug)]
pub struct MatplotlibFigure {
    /// Figure title.
    pub title: String,
    /// X-axis label.
    pub xlabel: String,
    /// Y-axis label.
    pub ylabel: String,
    /// Figure width in inches.
    pub width: f64,
    /// Figure height in inches.
    pub height: f64,
    /// Data series.
    pub series: Vec<DataSeries>,
    /// Whether to show grid.
    pub show_grid: bool,
    /// Legend location.
    pub legend_loc: String,
    /// DPI.
    pub dpi: u32,
}

/// A single data series for a Matplotlib plot.
#[derive(Clone, Debug)]
pub struct DataSeries {
    /// Series label.
    pub label: String,
    /// X data.
    pub x: Vec<f64>,
    /// Y data.
    pub y: Vec<f64>,
    /// Line style (e.g., "-", "--", "o").
    pub linestyle: String,
    /// Colour string.
    pub colour: String,
    /// Line width.
    pub linewidth: f64,
}

impl DataSeries {
    /// Create a new [`DataSeries`].
    pub fn new(label: impl Into<String>, x: Vec<f64>, y: Vec<f64>) -> Self {
        Self {
            label: label.into(),
            x,
            y,
            linestyle: "-".into(),
            colour: "blue".into(),
            linewidth: 1.5,
        }
    }
}

impl MatplotlibFigure {
    /// Create a new [`MatplotlibFigure`].
    pub fn new(
        title: impl Into<String>,
        xlabel: impl Into<String>,
        ylabel: impl Into<String>,
    ) -> Self {
        Self {
            title: title.into(),
            xlabel: xlabel.into(),
            ylabel: ylabel.into(),
            width: 8.0,
            height: 6.0,
            series: Vec::new(),
            show_grid: true,
            legend_loc: "best".into(),
            dpi: 150,
        }
    }

    /// Add a data series.
    pub fn add_series(&mut self, s: DataSeries) {
        self.series.push(s);
    }

    /// Export as a JSON descriptor string.
    pub fn to_json(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "title": "{}","#, self.title);
        let _ = writeln!(j, r#"  "xlabel": "{}","#, self.xlabel);
        let _ = writeln!(j, r#"  "ylabel": "{}","#, self.ylabel);
        let _ = writeln!(j, r#"  "figsize": [{}, {}],"#, self.width, self.height);
        let _ = writeln!(j, r#"  "dpi": {},"#, self.dpi);
        let _ = writeln!(j, r#"  "grid": {},"#, self.show_grid);
        let _ = writeln!(j, r#"  "legend_loc": "{}","#, self.legend_loc);
        let _ = writeln!(j, r#"  "series": ["#);
        for (i, s) in self.series.iter().enumerate() {
            let comma = if i + 1 < self.series.len() { "," } else { "" };
            let _ = writeln!(j, "    {{");
            let _ = writeln!(j, r#"      "label": "{}","#, s.label);
            let _ = writeln!(j, r#"      "linestyle": "{}","#, s.linestyle);
            let _ = writeln!(j, r#"      "color": "{}","#, s.colour);
            let _ = writeln!(j, r#"      "linewidth": {},"#, s.linewidth);
            let x_str: Vec<String> = s.x.iter().map(|v| v.to_string()).collect();
            let y_str: Vec<String> = s.y.iter().map(|v| v.to_string()).collect();
            let _ = writeln!(j, r#"      "x": [{}],"#, x_str.join(", "));
            let _ = writeln!(j, r#"      "y": [{}]"#, y_str.join(", "));
            let _ = writeln!(j, "    }}{}", comma);
        }
        let _ = writeln!(j, "  ]");
        let _ = writeln!(j, "}}");
        j
    }

    /// Generate a Python script to reproduce this figure.
    pub fn to_python_script(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "import matplotlib.pyplot as plt");
        let _ = writeln!(s, "import numpy as np");
        let _ = writeln!(
            s,
            "fig, ax = plt.subplots(figsize=({}, {}))",
            self.width, self.height
        );
        for series in &self.series {
            let x_vals: Vec<String> = series.x.iter().map(|v| v.to_string()).collect();
            let y_vals: Vec<String> = series.y.iter().map(|v| v.to_string()).collect();
            let _ = writeln!(
                s,
                "ax.plot([{}], [{}], '{}', color='{}', linewidth={}, label='{}')",
                x_vals.join(", "),
                y_vals.join(", "),
                series.linestyle,
                series.colour,
                series.linewidth,
                series.label
            );
        }
        let _ = writeln!(s, "ax.set_xlabel('{}')", self.xlabel);
        let _ = writeln!(s, "ax.set_ylabel('{}')", self.ylabel);
        let _ = writeln!(s, "ax.set_title('{}')", self.title);
        if self.show_grid {
            let _ = writeln!(s, "ax.grid(True)");
        }
        let _ = writeln!(s, "ax.legend(loc='{}')", self.legend_loc);
        let _ = writeln!(s, "plt.tight_layout()");
        let _ = writeln!(s, "plt.savefig('figure.png', dpi={})", self.dpi);
        s
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// D3.js data export
// ─────────────────────────────────────────────────────────────────────────────

/// A node in a force-directed graph.
#[derive(Clone, Debug)]
pub struct D3Node {
    /// Node ID.
    pub id: String,
    /// Node group (for colouring).
    pub group: u32,
    /// Optional radius.
    pub radius: f64,
    /// Optional label.
    pub label: String,
}

/// An edge in a force-directed graph.
#[derive(Clone, Debug)]
pub struct D3Edge {
    /// Source node ID.
    pub source: String,
    /// Target node ID.
    pub target: String,
    /// Edge value/weight.
    pub value: f64,
}

/// D3.js force-directed graph export.
#[derive(Clone, Debug, Default)]
pub struct D3ForceGraph {
    /// Graph nodes.
    pub nodes: Vec<D3Node>,
    /// Graph edges.
    pub edges: Vec<D3Edge>,
}

impl D3ForceGraph {
    /// Create a new [`D3ForceGraph`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node.
    pub fn add_node(&mut self, id: impl Into<String>, group: u32, radius: f64) {
        self.nodes.push(D3Node {
            id: id.into(),
            group,
            radius,
            label: String::new(),
        });
    }

    /// Add an edge.
    pub fn add_edge(&mut self, source: impl Into<String>, target: impl Into<String>, value: f64) {
        self.edges.push(D3Edge {
            source: source.into(),
            target: target.into(),
            value,
        });
    }

    /// Export to D3-compatible JSON.
    pub fn to_json(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "nodes": ["#);
        for (i, n) in self.nodes.iter().enumerate() {
            let comma = if i + 1 < self.nodes.len() { "," } else { "" };
            let _ = writeln!(
                j,
                r#"    {{"id": "{}", "group": {}, "radius": {}}}{}"#,
                n.id, n.group, n.radius, comma
            );
        }
        let _ = writeln!(j, r#"  ],"#);
        let _ = writeln!(j, r#"  "links": ["#);
        for (i, e) in self.edges.iter().enumerate() {
            let comma = if i + 1 < self.edges.len() { "," } else { "" };
            let _ = writeln!(
                j,
                r#"    {{"source": "{}", "target": "{}", "value": {}}}{}"#,
                e.source, e.target, e.value, comma
            );
        }
        let _ = writeln!(j, r#"  ]"#);
        let _ = writeln!(j, "}}");
        j
    }
}

/// Contour data for D3 contour plots.
#[derive(Clone, Debug)]
pub struct D3ContourData {
    /// Grid width.
    pub width: usize,
    /// Grid height.
    pub height: usize,
    /// Flat row-major data.
    pub values: Vec<f64>,
    /// Contour thresholds.
    pub thresholds: Vec<f64>,
    /// X extent \[xmin, xmax\].
    pub x_extent: [f64; 2],
    /// Y extent \[ymin, ymax\].
    pub y_extent: [f64; 2],
}

impl D3ContourData {
    /// Create a new [`D3ContourData`] instance.
    pub fn new(width: usize, height: usize, x_extent: [f64; 2], y_extent: [f64; 2]) -> Self {
        Self {
            width,
            height,
            values: vec![0.0; width * height],
            thresholds: Vec::new(),
            x_extent,
            y_extent,
        }
    }

    /// Set value at (ix, iy).
    pub fn set(&mut self, ix: usize, iy: usize, val: f64) {
        self.values[iy * self.width + ix] = val;
    }

    /// Add a contour threshold.
    pub fn add_threshold(&mut self, t: f64) {
        self.thresholds.push(t);
    }

    /// Auto-generate *n* evenly-spaced thresholds between min and max.
    pub fn auto_thresholds(&mut self, n: usize) {
        let vmin = self.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let vmax = self
            .values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        self.thresholds = (0..n)
            .map(|i| vmin + (vmax - vmin) * i as f64 / (n - 1).max(1) as f64)
            .collect();
    }

    /// Export to D3 contour-compatible JSON.
    pub fn to_json(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "width": {},"#, self.width);
        let _ = writeln!(j, r#"  "height": {},"#, self.height);
        let _ = writeln!(
            j,
            r#"  "x_extent": [{}, {}],"#,
            self.x_extent[0], self.x_extent[1]
        );
        let _ = writeln!(
            j,
            r#"  "y_extent": [{}, {}],"#,
            self.y_extent[0], self.y_extent[1]
        );
        let thresh_str: Vec<String> = self.thresholds.iter().map(|v| v.to_string()).collect();
        let _ = writeln!(j, r#"  "thresholds": [{}],"#, thresh_str.join(", "));
        let val_str: Vec<String> = self.values.iter().map(|v| v.to_string()).collect();
        let _ = writeln!(j, r#"  "values": [{}]"#, val_str.join(", "));
        let _ = writeln!(j, "}}");
        j
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WebGL buffer export
// ─────────────────────────────────────────────────────────────────────────────

/// WebGL-ready interleaved vertex buffer descriptor.
#[derive(Clone, Debug)]
pub struct WebGlBuffer {
    /// Buffer name.
    pub name: String,
    /// Interleaved float data: \[x, y, z, nx, ny, nz, u, v, ...\].
    pub data: Vec<f32>,
    /// Index buffer.
    pub indices: Vec<u32>,
    /// Stride in bytes.
    pub stride: u32,
    /// Attribute offsets: (name → byte offset).
    pub attributes: Vec<(String, u32, u32)>, // (name, offset, component_count)
}

impl WebGlBuffer {
    /// Create a new [`WebGlBuffer`] from separate arrays (positions + normals + uvs).
    pub fn from_mesh(
        name: impl Into<String>,
        positions: &[[f32; 3]],
        normals: &[[f32; 3]],
        uvs: &[[f32; 2]],
        indices: Vec<u32>,
    ) -> Self {
        let mut data = Vec::with_capacity(positions.len() * 8);
        for i in 0..positions.len() {
            data.extend_from_slice(&positions[i]);
            if i < normals.len() {
                data.extend_from_slice(&normals[i]);
            } else {
                data.extend_from_slice(&[0.0, 1.0, 0.0]);
            }
            if i < uvs.len() {
                data.extend_from_slice(&uvs[i]);
            } else {
                data.extend_from_slice(&[0.0, 0.0]);
            }
        }
        let stride = 8 * 4; // 8 floats × 4 bytes
        let attributes = vec![
            ("position".to_string(), 0u32, 3u32),
            ("normal".to_string(), 12u32, 3u32),
            ("uv".to_string(), 24u32, 2u32),
        ];
        Self {
            name: name.into(),
            data,
            indices,
            stride,
            attributes,
        }
    }

    /// Export metadata as JSON.
    pub fn to_json_meta(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "name": "{}","#, self.name);
        let _ = writeln!(j, r#"  "vertex_count": {},"#, self.data.len() / 8);
        let _ = writeln!(j, r#"  "index_count": {},"#, self.indices.len());
        let _ = writeln!(j, r#"  "stride": {},"#, self.stride);
        let _ = writeln!(j, r#"  "attributes": ["#);
        for (i, (aname, offset, count)) in self.attributes.iter().enumerate() {
            let comma = if i + 1 < self.attributes.len() {
                ","
            } else {
                ""
            };
            let _ = writeln!(
                j,
                r#"    {{"name": "{}", "offset": {}, "components": {}}}{}"#,
                aname, offset, count, comma
            );
        }
        let _ = writeln!(j, "  ]");
        let _ = writeln!(j, "}}");
        j
    }

    /// Export to binary (little-endian f32 array).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.data.len() * 4);
        for &v in &self.data {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        bytes
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.data.len() / 8
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// glTF 2.0 physics annotations
// ─────────────────────────────────────────────────────────────────────────────

/// Physics body descriptor for glTF KHR_physics_rigid_bodies extension.
#[derive(Clone, Debug)]
pub struct GltfPhysicsBody {
    /// Node name in the glTF scene.
    pub node_name: String,
    /// Mass in kg.
    pub mass: f64,
    /// Linear velocity \[vx, vy, vz\].
    pub linear_velocity: [f64; 3],
    /// Angular velocity \[wx, wy, wz\].
    pub angular_velocity: [f64; 3],
    /// Is this a static (non-moving) body?
    pub is_static: bool,
    /// Friction coefficient.
    pub friction: f64,
    /// Restitution (bounciness).
    pub restitution: f64,
    /// Collider shape type.
    pub collider: GltfColliderShape,
}

/// Collider shape for glTF physics.
#[derive(Clone, Debug)]
pub enum GltfColliderShape {
    /// Sphere collider with radius.
    Sphere(f64),
    /// Box collider with half-extents.
    Box([f64; 3]),
    /// Capsule with radius and half-height.
    Capsule(f64, f64),
    /// Convex hull (uses mesh).
    ConvexHull,
    /// Trimesh (uses mesh).
    TriMesh,
}

impl GltfPhysicsBody {
    /// Create a default dynamic sphere body.
    pub fn sphere(node_name: impl Into<String>, mass: f64, radius: f64) -> Self {
        Self {
            node_name: node_name.into(),
            mass,
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_static: false,
            friction: 0.5,
            restitution: 0.3,
            collider: GltfColliderShape::Sphere(radius),
        }
    }

    /// Serialize to glTF extension JSON fragment.
    pub fn to_gltf_json(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "KHR_physics_rigid_bodies": {{"#);
        let _ = writeln!(j, r#"    "motion": {{"#);
        let _ = writeln!(j, r#"      "isKinematic": {},"#, self.is_static);
        let _ = writeln!(j, r#"      "mass": {},"#, self.mass);
        let _ = writeln!(
            j,
            r#"      "linearVelocity": [{}, {}, {}],"#,
            self.linear_velocity[0], self.linear_velocity[1], self.linear_velocity[2]
        );
        let _ = writeln!(
            j,
            r#"      "angularVelocity": [{}, {}, {}]"#,
            self.angular_velocity[0], self.angular_velocity[1], self.angular_velocity[2]
        );
        let _ = writeln!(j, r#"    }},"#);
        let _ = writeln!(j, r#"    "collider": {{"#);
        match &self.collider {
            GltfColliderShape::Sphere(r) => {
                let _ = writeln!(j, r#"      "shape": "sphere","#);
                let _ = writeln!(j, r#"      "sphere": {{"radius": {}}}"#, r);
            }
            GltfColliderShape::Box(he) => {
                let _ = writeln!(j, r#"      "shape": "box","#);
                let _ = writeln!(
                    j,
                    r#"      "box": {{"halfExtents": [{}, {}, {}]}}"#,
                    he[0], he[1], he[2]
                );
            }
            GltfColliderShape::Capsule(r, h) => {
                let _ = writeln!(j, r#"      "shape": "capsule","#);
                let _ = writeln!(
                    j,
                    r#"      "capsule": {{"radius": {}, "height": {}}}"#,
                    r,
                    h * 2.0
                );
            }
            GltfColliderShape::ConvexHull => {
                let _ = writeln!(j, r#"      "shape": "convexHull""#);
            }
            GltfColliderShape::TriMesh => {
                let _ = writeln!(j, r#"      "shape": "trimesh""#);
            }
        }
        let _ = writeln!(j, r#"    }},"#);
        let _ = writeln!(
            j,
            r#"    "physicsMaterial": {{"friction": {}, "restitution": {}}}"#,
            self.friction, self.restitution
        );
        let _ = writeln!(j, r#"  }}"#);
        let _ = writeln!(j, "}}");
        j
    }
}

/// A complete glTF scene with physics annotations.
#[derive(Clone, Debug, Default)]
pub struct GltfPhysicsScene {
    /// Scene name.
    pub name: String,
    /// Physics bodies.
    pub bodies: Vec<GltfPhysicsBody>,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Fixed time step.
    pub fixed_dt: f64,
}

impl GltfPhysicsScene {
    /// Create a new [`GltfPhysicsScene`].
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            bodies: Vec::new(),
            gravity: [0.0, -9.81, 0.0],
            fixed_dt: 1.0 / 60.0,
        }
    }

    /// Add a physics body.
    pub fn add_body(&mut self, body: GltfPhysicsBody) {
        self.bodies.push(body);
    }

    /// Serialize scene physics extension to JSON.
    pub fn to_scene_json(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "name": "{}","#, self.name);
        let _ = writeln!(
            j,
            r#"  "gravity": [{}, {}, {}],"#,
            self.gravity[0], self.gravity[1], self.gravity[2]
        );
        let _ = writeln!(j, r#"  "fixedTimestep": {},"#, self.fixed_dt);
        let _ = writeln!(j, r#"  "bodies": ["#);
        for (i, body) in self.bodies.iter().enumerate() {
            let comma = if i + 1 < self.bodies.len() { "," } else { "" };
            let _ = write!(j, r#"    {{"node": "{}"}}{}"#, body.node_name, comma);
            let _ = writeln!(j);
        }
        let _ = writeln!(j, "  ]");
        let _ = writeln!(j, "}}");
        j
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VDB sparse volume export
// ─────────────────────────────────────────────────────────────────────────────

/// A sparse voxel in a VDB-like grid.
#[derive(Clone, Debug)]
pub struct VdbVoxel {
    /// Voxel grid coordinates.
    pub ijk: [i32; 3],
    /// Scalar value.
    pub value: f64,
}

/// A sparse VDB volume grid (minimal ASCII representation).
#[derive(Clone, Debug)]
pub struct VdbSparseGrid {
    /// Grid name.
    pub name: String,
    /// Grid type ("float", "vec3s", etc.).
    pub grid_type: String,
    /// World-space voxel size.
    pub voxel_size: f64,
    /// Background (default) value.
    pub background: f64,
    /// Populated voxels.
    pub voxels: Vec<VdbVoxel>,
}

impl VdbSparseGrid {
    /// Create a new [`VdbSparseGrid`].
    pub fn new(name: impl Into<String>, voxel_size: f64, background: f64) -> Self {
        Self {
            name: name.into(),
            grid_type: "float".into(),
            voxel_size,
            background,
            voxels: Vec::new(),
        }
    }

    /// Add a voxel.
    pub fn add_voxel(&mut self, i: i32, j: i32, k: i32, value: f64) {
        self.voxels.push(VdbVoxel {
            ijk: [i, j, k],
            value,
        });
    }

    /// Populate from a ScalarField3D (only non-background voxels).
    pub fn from_scalar_field(
        field: &ScalarField3D,
        threshold: f64,
        name: impl Into<String>,
    ) -> Self {
        let mut grid = Self::new(name, field.dx, 0.0);
        for ix in 0..field.nx {
            for iy in 0..field.ny {
                for iz in 0..field.nz {
                    let v = field.get(ix, iy, iz);
                    if v.abs() > threshold {
                        grid.add_voxel(ix as i32, iy as i32, iz as i32, v);
                    }
                }
            }
        }
        grid
    }

    /// Write a minimal ASCII VDB-like header.
    pub fn write_ascii_header(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "#VDB ASCII v1.0");
        let _ = writeln!(s, "name: {}", self.name);
        let _ = writeln!(s, "type: {}", self.grid_type);
        let _ = writeln!(s, "voxel_size: {}", self.voxel_size);
        let _ = writeln!(s, "background: {}", self.background);
        let _ = writeln!(s, "n_active: {}", self.voxels.len());
        for v in &self.voxels {
            let _ = writeln!(s, "v {} {} {} {}", v.ijk[0], v.ijk[1], v.ijk[2], v.value);
        }
        s
    }

    /// Number of active voxels.
    pub fn n_active(&self) -> usize {
        self.voxels.len()
    }

    /// Bounding box as \[(imin,imax),(jmin,jmax),(kmin,kmax)\].
    pub fn bbox(&self) -> [(i32, i32); 3] {
        if self.voxels.is_empty() {
            return [(0, 0); 3];
        }
        let mut lo = self.voxels[0].ijk;
        let mut hi = self.voxels[0].ijk;
        for v in &self.voxels {
            for d in 0..3 {
                if v.ijk[d] < lo[d] {
                    lo[d] = v.ijk[d];
                }
                if v.ijk[d] > hi[d] {
                    hi[d] = v.ijk[d];
                }
            }
        }
        [(lo[0], hi[0]), (lo[1], hi[1]), (lo[2], hi[2])]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OpenEXR multi-channel image descriptor
// ─────────────────────────────────────────────────────────────────────────────

/// A single channel in an OpenEXR image.
#[derive(Clone, Debug)]
pub struct ExrChannel {
    /// Channel name (e.g., "R", "G", "B", "depth", "normal.X").
    pub name: String,
    /// Channel data (half-float stored as f32).
    pub data: Vec<f32>,
    /// Pixel type ("FLOAT", "HALF", "UINT").
    pub pixel_type: String,
}

impl ExrChannel {
    /// Create a new [`ExrChannel`].
    pub fn new(name: impl Into<String>, width: usize, height: usize) -> Self {
        Self {
            name: name.into(),
            data: vec![0.0; width * height],
            pixel_type: "FLOAT".into(),
        }
    }

    /// Fill channel from a closure f(x, y) → f32.
    pub fn fill_from<F: Fn(usize, usize) -> f32>(&mut self, width: usize, height: usize, f: F) {
        for y in 0..height {
            for x in 0..width {
                self.data[y * width + x] = f(x, y);
            }
        }
    }
}

/// A multi-channel OpenEXR image descriptor.
#[derive(Clone, Debug)]
pub struct OpenExrImage {
    /// Image width.
    pub width: usize,
    /// Image height.
    pub height: usize,
    /// Channels.
    pub channels: Vec<ExrChannel>,
    /// Compression type.
    pub compression: String,
    /// Metadata key-value pairs.
    pub metadata: HashMap<String, String>,
}

impl OpenExrImage {
    /// Create a new [`OpenExrImage`].
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            channels: Vec::new(),
            compression: "PIZ".into(),
            metadata: HashMap::new(),
        }
    }

    /// Add a channel.
    pub fn add_channel(&mut self, ch: ExrChannel) {
        self.channels.push(ch);
    }

    /// Insert metadata.
    pub fn set_meta(&mut self, key: impl Into<String>, val: impl Into<String>) {
        self.metadata.insert(key.into(), val.into());
    }

    /// Write ASCII header descriptor.
    pub fn write_header(&self) -> String {
        let mut h = String::new();
        let _ = writeln!(h, "OPENEXR ASCII HEADER");
        let _ = writeln!(h, "width: {}", self.width);
        let _ = writeln!(h, "height: {}", self.height);
        let _ = writeln!(h, "compression: {}", self.compression);
        let _ = writeln!(h, "channels:");
        for ch in &self.channels {
            let _ = writeln!(h, "  {} ({})", ch.name, ch.pixel_type);
        }
        for (k, v) in &self.metadata {
            let _ = writeln!(h, "meta {} = {}", k, v);
        }
        h
    }

    /// Compute total bytes (f32) for all channels.
    pub fn total_bytes(&self) -> usize {
        self.channels.iter().map(|ch| ch.data.len() * 4).sum()
    }

    /// Create a standard RGBA image.
    pub fn rgba(width: usize, height: usize) -> Self {
        let mut img = Self::new(width, height);
        img.add_channel(ExrChannel::new("R", width, height));
        img.add_channel(ExrChannel::new("G", width, height));
        img.add_channel(ExrChannel::new("B", width, height));
        img.add_channel(ExrChannel::new("A", width, height));
        img
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cinema database format
// ─────────────────────────────────────────────────────────────────────────────

/// A Cinema Science Database (CDB) parameter definition.
#[derive(Clone, Debug)]
pub struct CinemaParameter {
    /// Parameter name.
    pub name: String,
    /// List of discrete values (stored as strings).
    pub values: Vec<String>,
    /// Whether this parameter is numeric.
    pub is_numeric: bool,
}

impl CinemaParameter {
    /// Create a new numeric parameter from a range.
    pub fn numeric_range(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            values: values.iter().map(|v| v.to_string()).collect(),
            is_numeric: true,
        }
    }

    /// Create a new string-valued parameter.
    pub fn string_param(name: impl Into<String>, values: Vec<String>) -> Self {
        Self {
            name: name.into(),
            values,
            is_numeric: false,
        }
    }
}

/// A Cinema Science Database (Spec D).
#[derive(Clone, Debug)]
pub struct CinemaDatabase {
    /// Database name.
    pub name: String,
    /// Parameters sweeping the parameter space.
    pub parameters: Vec<CinemaParameter>,
    /// Image file extension.
    pub file_ext: String,
    /// Entries: each entry is a map param_name → value + filename.
    pub entries: Vec<HashMap<String, String>>,
}

impl CinemaDatabase {
    /// Create a new [`CinemaDatabase`].
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            parameters: Vec::new(),
            file_ext: "png".into(),
            entries: Vec::new(),
        }
    }

    /// Add a parameter.
    pub fn add_parameter(&mut self, param: CinemaParameter) {
        self.parameters.push(param);
    }

    /// Add an image entry.
    pub fn add_entry(&mut self, param_values: HashMap<String, String>) {
        self.entries.push(param_values);
    }

    /// Compute total number of images (product of parameter cardinalities).
    pub fn total_images(&self) -> usize {
        self.parameters
            .iter()
            .map(|p| p.values.len().max(1))
            .product()
    }

    /// Write Cinema Spec D data.csv header.
    pub fn write_csv_header(&self) -> String {
        let mut cols: Vec<String> = self.parameters.iter().map(|p| p.name.clone()).collect();
        cols.push("FILE".into());
        cols.join(",")
    }

    /// Write Cinema Spec D data.csv body.
    pub fn write_csv_body(&self) -> String {
        let mut rows = Vec::new();
        for entry in &self.entries {
            let mut row_vals: Vec<String> = self
                .parameters
                .iter()
                .map(|p| entry.get(&p.name).cloned().unwrap_or_default())
                .collect();
            row_vals.push(entry.get("FILE").cloned().unwrap_or_default());
            rows.push(row_vals.join(","));
        }
        rows.join("\n")
    }

    /// Write full CSV (header + body).
    pub fn write_csv(&self) -> String {
        format!("{}\n{}", self.write_csv_header(), self.write_csv_body())
    }

    /// Generate a Cinema Spec D info.json.
    pub fn write_info_json(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "name_pattern": "{{FILE}}","#);
        let _ = writeln!(j, r#"  "arguments": {{"#);
        for (i, p) in self.parameters.iter().enumerate() {
            let comma = if i + 1 < self.parameters.len() {
                ","
            } else {
                ""
            };
            let vals: Vec<String> = p.values.iter().map(|v| format!(r#""{}""#, v)).collect();
            let _ = writeln!(j, r#"    "{}": [{}]{}"#, p.name, vals.join(", "), comma);
        }
        let _ = writeln!(j, r#"  }}"#);
        let _ = writeln!(j, "}}");
        j
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: colour map operations
// ─────────────────────────────────────────────────────────────────────────────

/// Built-in colour maps.
#[derive(Clone, Debug, PartialEq)]
pub enum ColourMap {
    /// Viridis.
    Viridis,
    /// Plasma.
    Plasma,
    /// Cool-warm (Bluered).
    CoolWarm,
    /// Greyscale.
    Greyscale,
    /// Hot (black-red-yellow-white).
    Hot,
}

impl ColourMap {
    /// Map scalar t ∈ \[0,1\] to RGBA \[r,g,b,a\] using this colour map.
    pub fn map(&self, t: f64) -> [f64; 4] {
        let t = t.clamp(0.0, 1.0);
        match self {
            ColourMap::Viridis => {
                // Approximate Viridis using polynomial coefficients
                let r = 0.267004 + t * (0.004874 + t * (0.329415 + t * (-0.001674)));
                let g = 0.004874 + t * (0.872325 + t * (-0.301631 + t * (-0.1)));
                let b = 0.329415 + t * (-0.635877 + t * (0.914499 + t * (-0.401658)));
                [r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0), 1.0]
            }
            ColourMap::Plasma => {
                let r = 0.050383 + t * (2.566707 + t * (-2.237019 + t * (0.816839)));
                let g = 0.029803 + t * (-0.390895 + t * (1.607541 + t * (-0.892576)));
                let b = 0.527975 + t * (1.016567 + t * (-2.476404 + t * (1.476605)));
                [r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0), 1.0]
            }
            ColourMap::CoolWarm => {
                let r = if t < 0.5 { 2.0 * t } else { 1.0 };
                let b = if t > 0.5 { 2.0 * (1.0 - t) } else { 1.0 };
                let g = 1.0 - 2.0 * (t - 0.5).abs();
                [r, g.max(0.0), b, 1.0]
            }
            ColourMap::Greyscale => [t, t, t, 1.0],
            ColourMap::Hot => {
                let r = (t * 3.0).min(1.0);
                let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
                let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
                [r, g, b, 1.0]
            }
        }
    }

    /// Name of the colour map.
    pub fn name(&self) -> &str {
        match self {
            ColourMap::Viridis => "viridis",
            ColourMap::Plasma => "plasma",
            ColourMap::CoolWarm => "coolwarm",
            ColourMap::Greyscale => "greyscale",
            ColourMap::Hot => "hot",
        }
    }
}

/// Convert a 2D scalar field slice to raw RGBA bytes using a colour map.
pub fn scalar_field_to_rgba(
    slice: &[f64],
    width: usize,
    height: usize,
    vmin: f64,
    vmax: f64,
    cmap: &ColourMap,
) -> Vec<u8> {
    let range = (vmax - vmin).max(1e-30);
    let mut bytes = Vec::with_capacity(width * height * 4);
    for &v in slice {
        let t = (v - vmin) / range;
        let rgba = cmap.map(t);
        for c in &rgba {
            bytes.push((c * 255.0).clamp(0.0, 255.0) as u8);
        }
    }
    bytes
}

// ─────────────────────────────────────────────────────────────────────────────
// Animation frame utilities
// ─────────────────────────────────────────────────────────────────────────────

/// A single animation frame carrying field data.
#[derive(Clone, Debug)]
pub struct AnimFrame {
    /// Time stamp.
    pub time: f64,
    /// Frame index.
    pub index: usize,
    /// Named scalar fields (name → data).
    pub scalars: HashMap<String, Vec<f64>>,
    /// Named vector fields (name → flat \[x,y,z, ...\]).
    pub vectors: HashMap<String, Vec<f64>>,
}

impl AnimFrame {
    /// Create a new [`AnimFrame`].
    pub fn new(time: f64, index: usize) -> Self {
        Self {
            time,
            index,
            scalars: HashMap::new(),
            vectors: HashMap::new(),
        }
    }

    /// Add a scalar field to this frame.
    pub fn add_scalar(&mut self, name: impl Into<String>, data: Vec<f64>) {
        self.scalars.insert(name.into(), data);
    }

    /// Add a vector field to this frame.
    pub fn add_vector(&mut self, name: impl Into<String>, data: Vec<f64>) {
        self.vectors.insert(name.into(), data);
    }
}

/// A sequence of animation frames.
#[derive(Clone, Debug, Default)]
pub struct AnimSequence {
    /// All frames.
    pub frames: Vec<AnimFrame>,
    /// Frame rate.
    pub fps: f64,
}

impl AnimSequence {
    /// Create a new [`AnimSequence`].
    pub fn new(fps: f64) -> Self {
        Self {
            frames: Vec::new(),
            fps,
        }
    }

    /// Push a frame.
    pub fn push(&mut self, frame: AnimFrame) {
        self.frames.push(frame);
    }

    /// Total duration in seconds.
    pub fn duration(&self) -> f64 {
        if self.frames.is_empty() {
            0.0
        } else {
            self.frames
                .last()
                .expect("collection should not be empty")
                .time
        }
    }

    /// Get frame at index.
    pub fn get_frame(&self, idx: usize) -> Option<&AnimFrame> {
        self.frames.get(idx)
    }

    /// Write a manifest JSON for the animation.
    pub fn write_manifest(&self) -> String {
        let mut j = String::new();
        let _ = writeln!(j, "{{");
        let _ = writeln!(j, r#"  "fps": {},"#, self.fps);
        let _ = writeln!(j, r#"  "n_frames": {},"#, self.frames.len());
        let _ = writeln!(j, r#"  "duration": {},"#, self.duration());
        if !self.frames.is_empty() {
            let scalar_names: Vec<String> = self.frames[0].scalars.keys().cloned().collect();
            let names_str: Vec<String> =
                scalar_names.iter().map(|n| format!(r#""{}""#, n)).collect();
            let _ = writeln!(j, r#"  "scalar_fields": [{}]"#, names_str.join(", "));
        }
        let _ = writeln!(j, "}}");
        j
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ScalarField3D ─────────────────────────────────────────────────────────

    #[test]
    fn test_scalar_field_set_get() {
        let mut f = ScalarField3D::new(4, 4, 4, 0.1, [0.0; 3]);
        f.set(1, 2, 3, 3.125);
        assert!((f.get(1, 2, 3) - 3.125).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_field_min_max() {
        let mut f = ScalarField3D::new(2, 2, 2, 0.1, [0.0; 3]);
        f.set(0, 0, 0, -5.0);
        f.set(1, 1, 1, 10.0);
        assert!((f.min_val() - (-5.0)).abs() < 1e-10);
        assert!((f.max_val() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_field_default_zero() {
        let f = ScalarField3D::new(3, 3, 3, 0.1, [0.0; 3]);
        assert_eq!(f.get(1, 1, 1), 0.0);
    }

    // ── ParaView state ────────────────────────────────────────────────────────

    #[test]
    fn test_paraview_state_contains_filename() {
        let cfg = ParaviewStateConfig::default_config("test.vtu");
        let s = write_paraview_state(&cfg);
        assert!(s.contains("test.vtu"));
    }

    #[test]
    fn test_paraview_state_contains_xml_header() {
        let cfg = ParaviewStateConfig::default_config("test.vtu");
        let s = write_paraview_state(&cfg);
        assert!(s.contains(r#"<?xml version="1.0"?>"#));
    }

    #[test]
    fn test_paraview_state_contains_camera() {
        let mut cfg = ParaviewStateConfig::default_config("test.vtu");
        cfg.camera_pos = [1.0, 2.0, 3.0];
        let s = write_paraview_state(&cfg);
        assert!(s.contains("CameraPosition"));
    }

    // ── VisIt database ────────────────────────────────────────────────────────

    #[test]
    fn test_visit_index_contains_files() {
        let mut db = VisItDatabase::new("sim", "vtk", 3);
        db.add_time(0.0);
        db.add_time(0.1);
        let idx = db.write_visit_index();
        assert!(idx.contains("sim000000.vtk"));
        assert!(idx.contains("sim000001.vtk"));
    }

    #[test]
    fn test_visit_manifest_json() {
        let mut db = VisItDatabase::new("run", "vtu", 3);
        db.add_time(0.0);
        db.add_variable("pressure");
        let m = db.write_manifest();
        assert!(m.contains("\"pressure\""));
        assert!(m.contains("\"n_times\": 1"));
    }

    // ── Blender export ────────────────────────────────────────────────────────

    #[test]
    fn test_blender_script_imports_bpy() {
        let tmpdir = std::env::temp_dir();
        let cfg = BlenderExportConfig::new("obj1", tmpdir.join("out.blend").to_str().unwrap_or(""));
        let s = cfg.generate_script("mesh.obj");
        assert!(s.contains("import bpy"));
    }

    #[test]
    fn test_blender_script_contains_object_name() {
        let tmpdir = std::env::temp_dir();
        let cfg =
            BlenderExportConfig::new("my_obj", tmpdir.join("out.blend").to_str().unwrap_or(""));
        let s = cfg.generate_script("mesh.obj");
        assert!(s.contains("my_obj"));
    }

    #[test]
    fn test_blender_volume_script() {
        let tmpdir = std::env::temp_dir();
        let cfg = BlenderExportConfig::new("vol", tmpdir.join("vol.blend").to_str().unwrap_or(""));
        let s = cfg.generate_volume_script("smoke.vdb");
        assert!(s.contains("smoke.vdb"));
        assert!(s.contains("ShaderNodeEmission"));
    }

    // ── Matplotlib JSON ───────────────────────────────────────────────────────

    #[test]
    fn test_matplotlib_json_contains_title() {
        let mut fig = MatplotlibFigure::new("My Plot", "X", "Y");
        fig.add_series(DataSeries::new("data", vec![1.0, 2.0], vec![3.0, 4.0]));
        let j = fig.to_json();
        assert!(j.contains("My Plot"));
    }

    #[test]
    fn test_matplotlib_python_script() {
        let fig = MatplotlibFigure::new("Test", "time", "value");
        let s = fig.to_python_script();
        assert!(s.contains("import matplotlib.pyplot as plt"));
    }

    #[test]
    fn test_matplotlib_series_data_in_json() {
        let mut fig = MatplotlibFigure::new("T", "x", "y");
        fig.add_series(DataSeries::new("series1", vec![0.0, 1.0], vec![2.0, 3.0]));
        let j = fig.to_json();
        assert!(j.contains("series1"));
    }

    // ── D3 force graph ────────────────────────────────────────────────────────

    #[test]
    fn test_d3_graph_json_nodes() {
        let mut g = D3ForceGraph::new();
        g.add_node("A", 1, 5.0);
        g.add_node("B", 2, 3.0);
        g.add_edge("A", "B", 1.0);
        let j = g.to_json();
        assert!(j.contains(r#""id": "A""#));
        assert!(j.contains(r#""source": "A""#));
    }

    #[test]
    fn test_d3_contour_json() {
        let mut cd = D3ContourData::new(4, 4, [0.0, 1.0], [0.0, 1.0]);
        cd.set(2, 2, 1.5);
        cd.auto_thresholds(5);
        let j = cd.to_json();
        assert!(j.contains("\"width\": 4"));
        assert!(!cd.thresholds.is_empty());
    }

    // ── WebGL buffer ──────────────────────────────────────────────────────────

    #[test]
    fn test_webgl_buffer_vertex_count() {
        let pos = vec![[0.0f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let nor = vec![[0.0f32, 0.0, 1.0]; 3];
        let uv = vec![[0.0f32, 0.0]; 3];
        let idx = vec![0u32, 1, 2];
        let buf = WebGlBuffer::from_mesh("tri", &pos, &nor, &uv, idx);
        assert_eq!(buf.vertex_count(), 3);
    }

    #[test]
    fn test_webgl_buffer_bytes_len() {
        let pos = vec![[0.0f32, 0.0, 0.0]];
        let nor = vec![[0.0f32, 0.0, 1.0]];
        let uv = vec![[0.0f32, 0.0]];
        let buf = WebGlBuffer::from_mesh("p", &pos, &nor, &uv, vec![0]);
        let bytes = buf.to_bytes();
        assert_eq!(bytes.len(), 8 * 4); // 8 floats × 4 bytes
    }

    #[test]
    fn test_webgl_buffer_meta_json() {
        let buf = WebGlBuffer::from_mesh("mesh", &[], &[], &[], vec![]);
        let j = buf.to_json_meta();
        assert!(j.contains("\"name\": \"mesh\""));
    }

    // ── glTF physics ──────────────────────────────────────────────────────────

    #[test]
    fn test_gltf_physics_body_json_sphere() {
        let body = GltfPhysicsBody::sphere("node1", 1.0, 0.5);
        let j = body.to_gltf_json();
        assert!(j.contains("sphere"));
        assert!(j.contains("KHR_physics_rigid_bodies"));
    }

    #[test]
    fn test_gltf_scene_json() {
        let mut scene = GltfPhysicsScene::new("test_scene");
        scene.add_body(GltfPhysicsBody::sphere("ball", 1.0, 0.5));
        let j = scene.to_scene_json();
        assert!(j.contains("test_scene"));
        assert!(j.contains("ball"));
    }

    // ── VDB sparse grid ───────────────────────────────────────────────────────

    #[test]
    fn test_vdb_add_voxel() {
        let mut g = VdbSparseGrid::new("smoke", 0.1, 0.0);
        g.add_voxel(1, 2, 3, 0.5);
        assert_eq!(g.n_active(), 1);
    }

    #[test]
    fn test_vdb_from_scalar_field() {
        let mut field = ScalarField3D::new(4, 4, 4, 0.1, [0.0; 3]);
        field.set(2, 2, 2, 5.0);
        let grid = VdbSparseGrid::from_scalar_field(&field, 1.0, "f");
        assert_eq!(grid.n_active(), 1);
    }

    #[test]
    fn test_vdb_ascii_header() {
        let mut g = VdbSparseGrid::new("density", 0.05, 0.0);
        g.add_voxel(0, 0, 0, 1.0);
        let h = g.write_ascii_header();
        assert!(h.contains("#VDB ASCII"));
        assert!(h.contains("density"));
    }

    #[test]
    fn test_vdb_bbox() {
        let mut g = VdbSparseGrid::new("vol", 0.1, 0.0);
        g.add_voxel(-1, 0, 0, 1.0);
        g.add_voxel(5, 3, 2, 2.0);
        let bb = g.bbox();
        assert_eq!(bb[0], (-1, 5));
    }

    // ── OpenEXR ───────────────────────────────────────────────────────────────

    #[test]
    fn test_exr_rgba_channels() {
        let img = OpenExrImage::rgba(8, 8);
        assert_eq!(img.channels.len(), 4);
    }

    #[test]
    fn test_exr_total_bytes() {
        let img = OpenExrImage::rgba(4, 4);
        assert_eq!(img.total_bytes(), 4 * 4 * 4 * 4); // 4 channels × 16 px × 4 bytes
    }

    #[test]
    fn test_exr_header_contains_channels() {
        let img = OpenExrImage::rgba(2, 2);
        let h = img.write_header();
        assert!(h.contains("R"));
        assert!(h.contains("G"));
    }

    #[test]
    fn test_exr_channel_fill() {
        let mut ch = ExrChannel::new("Z", 3, 3);
        ch.fill_from(3, 3, |x, y| (x + y) as f32);
        assert!((ch.data[4] - 2.0).abs() < 1e-6); // (1,1) → 2
    }

    // ── Cinema database ───────────────────────────────────────────────────────

    #[test]
    fn test_cinema_total_images() {
        let mut db = CinemaDatabase::new("sim");
        db.add_parameter(CinemaParameter::numeric_range(
            "phi",
            vec![0.0, 90.0, 180.0, 270.0],
        ));
        db.add_parameter(CinemaParameter::numeric_range(
            "theta",
            vec![0.0, 45.0, 90.0],
        ));
        assert_eq!(db.total_images(), 12);
    }

    #[test]
    fn test_cinema_csv_header() {
        let mut db = CinemaDatabase::new("test");
        db.add_parameter(CinemaParameter::numeric_range("time", vec![0.0, 1.0]));
        let h = db.write_csv_header();
        assert!(h.contains("time"));
        assert!(h.contains("FILE"));
    }

    #[test]
    fn test_cinema_info_json() {
        let mut db = CinemaDatabase::new("db");
        db.add_parameter(CinemaParameter::numeric_range("angle", vec![0.0, 90.0]));
        let j = db.write_info_json();
        assert!(j.contains("\"angle\""));
    }

    // ── Colour map ────────────────────────────────────────────────────────────

    #[test]
    fn test_colourmap_greyscale_zero() {
        let rgba = ColourMap::Greyscale.map(0.0);
        for &c in &rgba[..3] {
            assert!((c - 0.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_colourmap_greyscale_one() {
        let rgba = ColourMap::Greyscale.map(1.0);
        for &c in &rgba[..3] {
            assert!((c - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_colourmap_hot_red_at_quarter() {
        let rgba = ColourMap::Hot.map(0.4);
        assert!(rgba[0] > 0.5); // red
    }

    #[test]
    fn test_scalar_to_rgba_len() {
        let data = vec![0.0, 0.5, 1.0, 0.25];
        let bytes = scalar_field_to_rgba(&data, 2, 2, 0.0, 1.0, &ColourMap::CoolWarm);
        assert_eq!(bytes.len(), 16); // 4 pixels × 4 channels
    }

    // ── AnimSequence ──────────────────────────────────────────────────────────

    #[test]
    fn test_anim_sequence_duration() {
        let mut seq = AnimSequence::new(24.0);
        seq.push(AnimFrame::new(0.0, 0));
        seq.push(AnimFrame::new(1.0, 1));
        assert!((seq.duration() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_anim_frame_scalars() {
        let mut frame = AnimFrame::new(0.0, 0);
        frame.add_scalar("pressure", vec![1.0, 2.0, 3.0]);
        assert_eq!(frame.scalars["pressure"].len(), 3);
    }

    #[test]
    fn test_anim_manifest_json() {
        let mut seq = AnimSequence::new(30.0);
        let mut f = AnimFrame::new(0.0, 0);
        f.add_scalar("vel", vec![0.0]);
        seq.push(f);
        let m = seq.write_manifest();
        assert!(m.contains("\"n_frames\": 1"));
    }

    // ── Point3 ────────────────────────────────────────────────────────────────

    #[test]
    fn test_point3_dist() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(3.0, 4.0, 0.0);
        assert!((a.dist(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_point3_to_array() {
        let p = Point3::new(1.0, 2.0, 3.0);
        assert_eq!(p.to_array(), [1.0, 2.0, 3.0]);
    }
}
