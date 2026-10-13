// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Recast-equivalent navmesh construction from triangle soup.
//!
//! Pipeline: validate config → rasterize → filter walkable → compact →
//! erode by agent radius → regions → filter/merge regions → contours
//! (trace, simplify, bridge holes) → polymesh (ear clipping + convex merge) → NavMesh
//!
//! Units: `cell_size`, `cell_height`, `agent_*` and `edge_max_len` are world
//! units; `edge_max_error` is measured in cells; `region_min_size` and
//! `region_merge_size` are side lengths in cells (the area thresholds are their
//! squares, as in Recast).
//!
//! **v0.1 limitations:** no detail mesh, no tile partitioning, no off-mesh
//! connections, a single area type.  Corridors narrower than
//! `2 * agent_radius` are removed by the erosion stage by design.

pub mod compact;
pub mod contour;
pub mod erode;
pub(crate) mod grid;
pub mod polymesh;
pub mod rasterize;
pub mod region;
pub mod region_merge;
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
    /// Agent radius; the walkable area is eroded by `ceil(agent_radius / cell_size)`
    /// cells so that agents keep this clearance from walls. Default: 0.6
    pub agent_radius: f64,
    /// Maximum height step an agent can climb. Default: 0.9
    pub agent_max_climb: f64,
    /// Maximum slope angle (degrees) walkable surfaces can have. Default: 45.0
    pub agent_max_slope: f64,
    /// Minimum number of cells per region. Default: 8
    pub region_min_size: usize,
    /// Regions with at most `region_merge_size²` cells (or not touching any wall)
    /// are merged into a neighbouring region. Default: 20
    pub region_merge_size: usize,
    /// Maximum length of a contour wall edge in world units (0 = unlimited).
    /// Default: 12.0
    pub edge_max_len: f64,
    /// Maximum deviation (in cells) of the simplified contour from the raw
    /// region boundary. Default: 1.3
    pub edge_max_error: f64,
    /// Maximum vertices per polygon (3..=6). Default: 6
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

impl RecastConfig {
    /// Check that every parameter is usable by the pipeline.
    ///
    /// * `cell_size`, `cell_height` and `agent_height` must be finite and `> 0`;
    /// * `agent_radius`, `agent_max_climb`, `edge_max_len` and `edge_max_error`
    ///   must be finite and `>= 0`;
    /// * `agent_max_slope` must lie in `[0, 90)` degrees;
    /// * `max_verts_per_poly` must lie in `3..=6`;
    /// * the derived voxel counts (agent height / climb / radius) must fit the
    ///   `u16` voxel range and the region area thresholds must not overflow.
    pub fn validate(&self) -> Result<(), RecastError> {
        let invalid = |msg: String| Err(RecastError::InvalidConfig(msg));
        for (name, value) in [
            ("cell_size", self.cell_size),
            ("cell_height", self.cell_height),
            ("agent_height", self.agent_height),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return invalid(format!("{name} must be finite and > 0, got {value}"));
            }
        }
        for (name, value) in [
            ("agent_radius", self.agent_radius),
            ("agent_max_climb", self.agent_max_climb),
            ("edge_max_len", self.edge_max_len),
            ("edge_max_error", self.edge_max_error),
        ] {
            if !(value.is_finite() && value >= 0.0) {
                return invalid(format!("{name} must be finite and >= 0, got {value}"));
            }
        }
        if !(self.agent_max_slope.is_finite() && (0.0..90.0).contains(&self.agent_max_slope)) {
            return invalid(format!(
                "agent_max_slope must be in [0, 90) degrees, got {}",
                self.agent_max_slope
            ));
        }
        if !(3..=polymesh::MAX_VERTS_PER_POLY).contains(&self.max_verts_per_poly) {
            return invalid(format!(
                "max_verts_per_poly must be in 3..={}, got {}",
                polymesh::MAX_VERTS_PER_POLY,
                self.max_verts_per_poly
            ));
        }
        let voxel_limit = f64::from(u16::MAX);
        for (name, voxels) in [
            (
                "agent_height / cell_height",
                self.agent_height / self.cell_height,
            ),
            (
                "agent_max_climb / cell_height",
                self.agent_max_climb / self.cell_height,
            ),
            (
                "agent_radius / cell_size",
                self.agent_radius / self.cell_size,
            ),
        ] {
            if !(voxels.is_finite() && voxels.ceil() < voxel_limit) {
                return invalid(format!(
                    "{name} = {voxels} voxels exceeds the u16 voxel range"
                ));
            }
        }
        for (name, side) in [
            ("region_min_size", self.region_min_size),
            ("region_merge_size", self.region_merge_size),
        ] {
            if side.checked_mul(side).is_none() {
                return invalid(format!("{name}² overflows usize (got {side})"));
            }
        }
        Ok(())
    }
}

/// Error type for Recast pipeline.
#[derive(Debug)]
pub enum RecastError {
    EmptyInput,
    /// The [`RecastConfig`] contains an unusable parameter.
    InvalidConfig(String),
    RasterizationFailed(String),
    RegionBuildFailed(String),
    ContourBuildFailed(String),
    PolyMeshBuildFailed(String),
}

impl std::fmt::Display for RecastError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RecastError::EmptyInput => write!(f, "input triangle soup is empty"),
            RecastError::InvalidConfig(s) => write!(f, "invalid recast config: {s}"),
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
        self.config.validate()?;
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
        let mut chf = compact::build_compact_heightfield(&hf, &self.config);

        // Stage 3b: Erode the walkable area by the agent radius, then rebuild
        // the compact heightfield (connectivity + distance field) without the
        // eroded spans.
        let walkable_radius =
            erode::walkable_radius_cells(self.config.agent_radius, self.config.cell_size);
        if walkable_radius > 0 && erode::erode_walkable_area(&mut hf, &chf, walkable_radius) > 0 {
            chf = compact::build_compact_heightfield(&hf, &self.config);
        }

        // Stage 4: Build regions, then drop tiny regions and merge small ones
        let mut rhf =
            region::build_regions(&chf, &self.config).map_err(RecastError::RegionBuildFailed)?;
        let min_region_area = self
            .config
            .region_min_size
            .saturating_mul(self.config.region_min_size);
        let merge_region_area = self
            .config
            .region_merge_size
            .saturating_mul(self.config.region_merge_size);
        region_merge::merge_and_filter_regions(&mut rhf, min_region_area, merge_region_area)
            .map_err(RecastError::RegionBuildFailed)?;

        // Stage 5: Build contours
        let cs = contour::build_contours(&rhf, &chf, &self.config)
            .map_err(RecastError::ContourBuildFailed)?;

        // Stage 6: Build polygon mesh
        let pm = polymesh::build_poly_mesh(&cs, &self.config)
            .map_err(RecastError::PolyMeshBuildFailed)?;

        Ok(pm)
    }
}
