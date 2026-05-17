// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Colour helpers for the debug-tools surface.
//!
//! Re-exports the [`Rgba`] type used by every overlay primitive so JavaScript
//! can construct colour values without a heavy interop layer.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::to_js_value;

/// RGBA colour represented as four u8 components.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rgba {
    /// Red channel (0–255).
    pub r: u8,
    /// Green channel (0–255).
    pub g: u8,
    /// Blue channel (0–255).
    pub b: u8,
    /// Alpha channel (0–255).
    pub a: u8,
}

impl Rgba {
    /// Construct from components.
    pub fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }

    /// Opaque red.
    pub fn red() -> Self {
        Self::new(255, 0, 0, 255)
    }
    /// Opaque green.
    pub fn green() -> Self {
        Self::new(0, 255, 0, 255)
    }
    /// Opaque blue.
    pub fn blue() -> Self {
        Self::new(0, 0, 255, 255)
    }
    /// Opaque white.
    pub fn white() -> Self {
        Self::new(255, 255, 255, 255)
    }
    /// Opaque yellow.
    pub fn yellow() -> Self {
        Self::new(255, 255, 0, 255)
    }
    /// Transparent black.
    pub fn transparent() -> Self {
        Self::new(0, 0, 0, 0)
    }
}

#[wasm_bindgen]
impl Rgba {
    /// Construct from RGBA components (JS-compatible).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(r: u8, g: u8, b: u8, a: u8) -> Rgba {
        Rgba::new(r, g, b, a)
    }

    /// Opaque red preset.
    #[wasm_bindgen(js_name = "red")]
    pub fn red_js() -> Rgba {
        Rgba::red()
    }

    /// Opaque green preset.
    #[wasm_bindgen(js_name = "green")]
    pub fn green_js() -> Rgba {
        Rgba::green()
    }

    /// Opaque blue preset.
    #[wasm_bindgen(js_name = "blue")]
    pub fn blue_js() -> Rgba {
        Rgba::blue()
    }

    /// Opaque white preset.
    #[wasm_bindgen(js_name = "white")]
    pub fn white_js() -> Rgba {
        Rgba::white()
    }

    /// Opaque yellow preset.
    #[wasm_bindgen(js_name = "yellow")]
    pub fn yellow_js() -> Rgba {
        Rgba::yellow()
    }

    /// Fully transparent preset.
    #[wasm_bindgen(js_name = "transparent")]
    pub fn transparent_js() -> Rgba {
        Rgba::transparent()
    }

    /// Serialise to a `JsValue` object `{ r, g, b, a }`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
