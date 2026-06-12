// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Naga front-end validation for the crate's WGSL kernels.
//!
//! Unlike the substring checks in `shaders::functions::validate_wgsl_structure`,
//! these tests run the real Pure-Rust naga WGSL parser and validator over each
//! shader source. A shader that fails to parse or type-check (e.g. the legacy
//! `pub(super)` leakage that is invalid WGSL, or a malformed expression) makes
//! the corresponding test fail. No GPU adapter is required — this is a pure
//! compile-time validation of the WGSL itself.

/// Parse + validate one WGSL source through naga; panics with the shader name
/// and the naga diagnostic on failure.
fn validate_wgsl(name: &str, src: &str) {
    let module = match naga::front::wgsl::parse_str(src) {
        Ok(m) => m,
        Err(e) => panic!("WGSL parse failed for `{name}`:\n{}", e.emit_to_string(src)),
    };
    let mut validator = naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::all(),
    );
    if let Err(e) = validator.validate(&module) {
        panic!("WGSL validation failed for `{name}`: {e:?}");
    }
}

// ── SPH pipeline kernels (the dispatched deliverable) ─────────────────────────

#[test]
fn sph_cell_list_wgsl_compiles() {
    validate_wgsl(
        "SPH_CELL_LIST_WGSL",
        oxiphysics_gpu::kernels_wgsl::SPH_CELL_LIST_WGSL,
    );
}

#[test]
fn sph_density_grid_wgsl_compiles() {
    validate_wgsl(
        "SPH_DENSITY_GRID_WGSL",
        oxiphysics_gpu::kernels_wgsl::SPH_DENSITY_GRID_WGSL,
    );
}

#[test]
fn sph_force_grid_wgsl_compiles() {
    validate_wgsl(
        "SPH_FORCE_GRID_WGSL",
        oxiphysics_gpu::kernels_wgsl::SPH_FORCE_GRID_WGSL,
    );
}

#[test]
fn sph_integrate_wgsl_compiles() {
    validate_wgsl(
        "SPH_INTEGRATE_WGSL",
        oxiphysics_gpu::kernels_wgsl::SPH_INTEGRATE_WGSL,
    );
}

#[test]
fn sph_boundary_wgsl_compiles() {
    validate_wgsl(
        "SPH_BOUNDARY_WGSL",
        oxiphysics_gpu::kernels_wgsl::SPH_BOUNDARY_WGSL,
    );
}

// ── Reduction-primitive kernels (already dispatched elsewhere) ────────────────

#[test]
fn reduction_kernels_compile() {
    use oxiphysics_gpu::kernels_wgsl::*;
    validate_wgsl("SCAN_WGSL", SCAN_WGSL);
    validate_wgsl("SCAN_ADD_WGSL", SCAN_ADD_WGSL);
    validate_wgsl("REDUCE_WGSL", REDUCE_WGSL);
    validate_wgsl("COMPACT_WGSL", COMPACT_WGSL);
    validate_wgsl("HISTOGRAM_WGSL", HISTOGRAM_WGSL);
    validate_wgsl("RADIX_HISTOGRAM_WGSL", RADIX_HISTOGRAM_WGSL);
    validate_wgsl("RADIX_SCATTER_WGSL", RADIX_SCATTER_WGSL);
    validate_wgsl("RADIX_SCATTER_PAIRS_WGSL", RADIX_SCATTER_PAIRS_WGSL);
}

// ── shader_registry templates (the 4 promoted + the rest) ─────────────────────

/// All ten built-in `shaders::functions` WGSL templates must now be valid WGSL.
///
/// This is the regression guard for the `pub(super)` leakage that previously
/// made `CELL_LIST_WGSL`, `SPH_FORCE_WGSL`, `INTEGRATE_WGSL`, and
/// `BOUNDARY_ENFORCE_WGSL` (plus six siblings) parse-fail under any real WGSL
/// compiler while passing the old substring-only validator.
#[test]
fn shader_registry_templates_compile() {
    use oxiphysics_gpu::shaders::{
        BOUNDARY_ENFORCE_WGSL, BROADPHASE_SORT_SHADER, CELL_LIST_WGSL, INTEGRATE_WGSL,
        LBM_BGK_D2Q9_WGSL, LBM_STREAMING_SHADER, RIGID_INTEGRATE_SHADER, SDF_COMPUTE_WGSL,
        SPH_DENSITY_WGSL, SPH_FORCE_WGSL,
    };
    validate_wgsl("SPH_DENSITY_WGSL", SPH_DENSITY_WGSL);
    validate_wgsl("INTEGRATE_WGSL", INTEGRATE_WGSL);
    validate_wgsl("LBM_BGK_D2Q9_WGSL", LBM_BGK_D2Q9_WGSL);
    validate_wgsl("CELL_LIST_WGSL", CELL_LIST_WGSL);
    validate_wgsl("SDF_COMPUTE_WGSL", SDF_COMPUTE_WGSL);
    validate_wgsl("LBM_STREAMING_SHADER", LBM_STREAMING_SHADER);
    validate_wgsl("RIGID_INTEGRATE_SHADER", RIGID_INTEGRATE_SHADER);
    validate_wgsl("BROADPHASE_SORT_SHADER", BROADPHASE_SORT_SHADER);
    validate_wgsl("SPH_FORCE_WGSL", SPH_FORCE_WGSL);
    validate_wgsl("BOUNDARY_ENFORCE_WGSL", BOUNDARY_ENFORCE_WGSL);
}
