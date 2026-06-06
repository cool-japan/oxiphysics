//! Auto-generated module structure

pub mod batchnarrowphase_traits;
pub mod contactfilter_traits;
pub mod dispatch;
pub mod epa;
pub mod functions;
pub mod functions_2;
pub mod gjk;
pub mod gjk_cache;
pub mod manifold;
pub mod obb_sat;
pub mod primitives;
pub mod simple;
pub mod specialized;
pub mod types;

// Re-export from dispatch (primary for CompoundShape, ShapeType, generate_speculative_contacts)
pub use dispatch::*;

// Re-export from epa (primary for MAX_ITERATIONS, TOLERANCE)
pub use epa::*;

// Re-export from functions/functions_2
pub use functions::*;
pub use functions_2::*;

// Re-export from gjk (primary for gjk_distance, SupportCache, GjkDistanceResult, RayCastResult)
pub use gjk::*;

// Re-export from gjk_cache: only unique items (skip duplicates: ShapeType, SupportCache, GjkDistanceResult, gjk_distance)
pub use gjk_cache::{
    BatchProximityEntry, CachedSupport, CsoPoint, EvictionPolicy, GjkCache, GjkCacheRegistry,
    GjkCacheStats, GjkContactPair, GjkPairCache, GjkStats, GjkTermination, GjkWarmStart,
    HitRateStats, JohnsonSubResult, NarrowphaseResult, PositionedGjkCache, ProximityResult,
    ShapeDesc, SimplexCache, TimestampedCacheEntry, TimestampedGjkRegistry, WarmStartedGjk,
    cache_still_valid, dispatch_gjk_query, gjk_lower_bound, gjk_proximity, gjk_reduce_simplex,
    johnson_closest_point, simplex_vertex_max_dist, simplex_vertex_min_dist,
    sphere_sphere_narrowphase, support_minkowski_diff,
};

// Re-export from manifold: only unique items (skip: generate_speculative_contacts)
pub use manifold::{
    ContactManifold, ContactPoint, PersistentPatch, SpeculativeContactPoint, WarmContactPoint,
    WarmManifold, WarmStartData, build_4point_manifold, build_box_manifold,
    clip_polygon_by_halfspace, coulomb_friction_limit, decompose_velocity, friction_basis,
    manifold_quality, project_point_to_plane, reduce_manifold, rotate_manifold_points,
    smooth_contact_normal, smooth_manifold_normals, sutherland_hodgman_2d, transfer_warm_start,
    unproject_point_from_plane,
};

pub use obb_sat::*;
// NOTE: primitives::ClosestFeature conflicts with gjk::ClosestFeature.
// We skip the glob re-export for primitives; use narrowphase::primitives:: explicitly.

// Re-export from types: only unique items (skip: CompoundShape, RayCastResult)
pub use types::{
    BatchNarrowPhase, ContactFeature, ContactFilter, FeatureContact, NarrowPhaseContact,
    PointQueryResult, SegmentCastResult, ShapeKind, TriangleMesh,
};
