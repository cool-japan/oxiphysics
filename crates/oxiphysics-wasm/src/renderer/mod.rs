// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Canvas renderer hints for browser-based physics visualizations.
//!
//! This module provides lightweight data structures used to communicate
//! per-frame rendering information from the physics engine to a JavaScript
//! canvas 2D or WebGL renderer. All types are plain-data and serialize cleanly
//! to JSON for `postMessage` transfer or structured clone.
//!
//! ## Submodules
//!
//! - `core`: core hint types (AABB, contacts, velocity vectors, per-body hints)
//! - `wasm_renderer`: `WasmRenderer` with SSAO/instanced mesh / screen-space AABB
//! - `pipeline`: GPU pipeline descriptor helpers (blend modes, G-buffers, frame graph)

pub mod core;
pub mod pipeline;
pub mod wasm_renderer;

// Re-export everything so callers can continue to use `crate::renderer::Aabb`, etc.
pub use core::{
    Aabb, BodyRenderHint, ContactPointHint, DebugDrawCommand, RendererHints, VelocityVector,
    build_draw_commands,
};
pub use pipeline::{
    BlendMode, DepthCompare, FrameGraph, GBufferLayout, InstanceData, ParticleRenderPass,
    PostProcessPass, PostProcessPipeline, RenderPipelineDesc, ShadowMapPass,
};
pub use wasm_renderer::{InstancedMeshEntry, ScreenSpaceAabb, SsaoConfig, WasmRenderer};
