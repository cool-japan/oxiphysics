// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Post-processing effects: depth of field, bloom, tone mapping, and vignette.
//!
//! All types operate on a CPU-side framebuffer ([`Image`]) and use plain
//! `f32`/`f64` arrays — no external linear-algebra dependencies.

pub mod advanced_effects;
pub mod core;
pub mod effects;
pub mod filters;
pub mod image_filters;
pub mod processor;

// Re-export everything for backwards compatibility
pub use self::core::*;
pub use advanced_effects::*;
pub use effects::*;
pub use filters::*;
pub use image_filters::*;
pub use processor::*;
