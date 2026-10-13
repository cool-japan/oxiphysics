// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Watershed region partitioning of the compact heightfield.
//!
//! [`build_regions`] labels every walkable span with a region id.  The regions
//! are the basins of the distance field (see the `watershed` submodule):
//! seeds start at the local maxima of the distance to the nearest wall and
//! grow level by level, so a doorway or a pillar splits a room into several
//! regions.  The raw partition is usually finer than wanted;
//! [`merge_and_filter_regions`](super::region_merge::merge_and_filter_regions)
//! merges small regions into their neighbours afterwards.

use super::RecastConfig;
use super::compact::CompactHeightfield;

mod watershed;

/// Compact heightfield with region IDs assigned to each span.
pub struct RegionHeightfield {
    pub compact: CompactHeightfield,
}

/// Partition the walkable spans of `chf` into watershed regions.
///
/// Every span of the returned heightfield carries its region id in
/// `CompactSpan::reg` (`1..=n`, `0` = no region).  Connected groups of walkable
/// spans with fewer than `region_min_size²` spans can never form a region and
/// stay at `0`.
///
/// # Errors
///
/// Fails when the field needs more than 65535 regions (the ids are 16-bit) or
/// when the compact heightfield is inconsistent (distance field of the wrong
/// length, more spans than a 32-bit index can address).
pub fn build_regions(
    chf: &CompactHeightfield,
    cfg: &RecastConfig,
) -> Result<RegionHeightfield, String> {
    let min_region_area = cfg.region_min_size.saturating_mul(cfg.region_min_size);
    let regions = watershed::partition(chf, min_region_area)?;
    let mut compact = clone_compact_hf(chf);
    for (span, &region) in compact.spans.iter_mut().zip(&regions) {
        span.reg = region;
    }
    Ok(RegionHeightfield { compact })
}

fn clone_compact_hf(chf: &CompactHeightfield) -> CompactHeightfield {
    CompactHeightfield {
        width: chf.width,
        height: chf.height,
        span_count: chf.span_count,
        walkable_height: chf.walkable_height,
        walkable_climb: chf.walkable_climb,
        bmin: chf.bmin,
        bmax: chf.bmax,
        cs: chf.cs,
        ch: chf.ch,
        cells: chf.cells.clone(),
        spans: chf.spans.clone(),
        dist: chf.dist.clone(),
        max_dist: chf.max_dist,
    }
}
