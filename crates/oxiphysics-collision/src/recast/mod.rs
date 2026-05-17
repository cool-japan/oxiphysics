// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Recast-equivalent navmesh construction from triangle soup.
//!
//! Pipeline: rasterize → filter walkable → compact → regions → contours → polymesh → NavMesh
//!
//! **v0.1 limitations:** no detail mesh, no tile partitioning, no off-mesh connections.
//! Long narrow corridors (width < 2*cell_size) may produce empty regions (known limitation).

pub mod compact;
pub mod contour;
pub mod polymesh;
pub mod rasterize;
pub mod region;
pub mod walkable;

use crate::recast::polymesh::PolyMesh;

/// Configuration for the Recast pipeline.
#[derive(Debug, Clone)]
pub struct RecastConfig {
    /// Voxel size in world units (XZ dimensions). Default: 0.3
    pub cell_size: f64,
    /// Voxel height in world units (Y dimension). Default: 0.2
    pub cell_height: f64,
    /// Agent height (minimum clearance above walkable surface). Default: 2.0
    pub agent_height: f64,
    /// Agent radius. Default: 0.6
    pub agent_radius: f64,
    /// Maximum height step an agent can climb. Default: 0.9
    pub agent_max_climb: f64,
    /// Maximum slope angle (degrees) walkable surfaces can have. Default: 45.0
    pub agent_max_slope: f64,
    /// Minimum number of cells per region. Default: 8
    pub region_min_size: usize,
    /// Size threshold for region merging. Default: 20
    pub region_merge_size: usize,
    /// Maximum contour edge length. Default: 12.0
    pub edge_max_len: f64,
    /// Maximum contour simplification error. Default: 1.3
    pub edge_max_error: f64,
    /// Maximum vertices per polygon. Default: 6
    pub max_verts_per_poly: usize,
}

impl Default for RecastConfig {
    fn default() -> Self {
        Self {
            cell_size: 0.3,
            cell_height: 0.2,
            agent_height: 2.0,
            agent_radius: 0.6,
            agent_max_climb: 0.9,
            agent_max_slope: 45.0,
            region_min_size: 8,
            region_merge_size: 20,
            edge_max_len: 12.0,
            edge_max_error: 1.3,
            max_verts_per_poly: 6,
        }
    }
}

/// Error type for Recast pipeline.
#[derive(Debug)]
pub enum RecastError {
    EmptyInput,
    RasterizationFailed(String),
    RegionBuildFailed(String),
    ContourBuildFailed(String),
    PolyMeshBuildFailed(String),
}

impl std::fmt::Display for RecastError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RecastError::EmptyInput => write!(f, "input triangle soup is empty"),
            RecastError::RasterizationFailed(s) => write!(f, "rasterization failed: {s}"),
            RecastError::RegionBuildFailed(s) => write!(f, "region build failed: {s}"),
            RecastError::ContourBuildFailed(s) => write!(f, "contour build failed: {s}"),
            RecastError::PolyMeshBuildFailed(s) => write!(f, "poly mesh build failed: {s}"),
        }
    }
}

impl std::error::Error for RecastError {}

/// The main navmesh builder.
pub struct RecastBuilder {
    pub config: RecastConfig,
}

impl RecastBuilder {
    pub fn new(config: RecastConfig) -> Self {
        Self { config }
    }

    pub fn with_default_config() -> Self {
        Self::new(RecastConfig::default())
    }

    /// Run the full pipeline and return a `PolyMesh`.
    ///
    /// The caller is responsible for converting `PolyMesh` to the application's
    /// `NavMesh` type using `polymesh::poly_mesh_to_nav_mesh_primitives`.
    pub fn build(&self, tris: &[[f64; 9]]) -> Result<PolyMesh, RecastError> {
        if tris.is_empty() {
            return Err(RecastError::EmptyInput);
        }

        // Stage 1: Rasterize into heightfield
        let mut hf = rasterize::rasterize_triangles(tris, &self.config)
            .map_err(RecastError::RasterizationFailed)?;

        // Stage 2: Filter walkable spans
        walkable::filter_low_hanging_obstacles(&mut hf, &self.config);
        walkable::filter_ledge_spans(&mut hf, &self.config);
        walkable::filter_walkable_low_height(&mut hf, &self.config);

        // Stage 3: Compact heightfield
        let chf = compact::build_compact_heightfield(&hf, &self.config);

        // Stage 4: Build regions
        let rhf =
            region::build_regions(&chf, &self.config).map_err(RecastError::RegionBuildFailed)?;

        // Stage 5: Build contours
        let cs = contour::build_contours(&rhf, &chf, &self.config)
            .map_err(RecastError::ContourBuildFailed)?;

        // Stage 6: Build polygon mesh
        let pm = polymesh::build_poly_mesh(&cs, &self.config)
            .map_err(RecastError::PolyMeshBuildFailed)?;

        Ok(pm)
    }
}
