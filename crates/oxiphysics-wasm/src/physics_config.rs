// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics configuration types for browser-based WASM simulations.
//!
//! This module provides structured configuration for controlling gravity,
//! substep counts, broadphase type, solver parameters, and other global
//! simulation settings. All types serialize to JSON for storage or transfer.
//!
//! ## Typical usage
//!
//! ```no_run
//! use oxiphysics_wasm::physics_config::{PhysicsConfig, BroadphaseType, GravityPreset};
//!
//! let cfg = PhysicsConfig::builder()
//!     .gravity_preset(GravityPreset::Earth)
//!     .broadphase(BroadphaseType::SweepAndPrune)
//!     .substeps(4)
//!     .build();
//!
//! assert_eq!(cfg.substeps, 4);
//! ```

use crate::types::SimulationConfig;
use crate::wasm_helpers::{err_to_jsvalue, to_js_value};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// BroadphaseType
// ---------------------------------------------------------------------------

/// Broadphase algorithm selector.
///
/// The broadphase reduces the number of narrowphase collision tests by
/// quickly identifying pairs of bodies that *might* be colliding.
///
/// All variants are unit-only so this enum can be exposed directly to JS.
#[wasm_bindgen(js_name = "PhysicsBroadphaseType")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum BroadphaseType {
    /// Brute-force O(n²) pair testing. Fine for small scenes (<50 bodies).
    BruteForce,
    /// Sort-and-sweep along the X axis. Good for sparse, axis-aligned scenes.
    #[default]
    SweepAndPrune,
    /// Uniform spatial grid. Efficient when bodies are roughly evenly distributed.
    UniformGrid,
    /// BVH (bounding volume hierarchy). Best for large, dynamic scenes.
    Bvh,
    /// Dynamic AABB tree (like Bullet/Box2D). Incremental updates.
    DynamicAabbTree,
}

impl BroadphaseType {
    /// Return a human-readable name for this broadphase type.
    pub fn name(&self) -> &'static str {
        match self {
            Self::BruteForce => "BruteForce",
            Self::SweepAndPrune => "SweepAndPrune",
            Self::UniformGrid => "UniformGrid",
            Self::Bvh => "Bvh",
            Self::DynamicAabbTree => "DynamicAabbTree",
        }
    }

    /// Returns `true` if this broadphase supports dynamic body updates efficiently.
    pub fn supports_dynamic_updates(&self) -> bool {
        matches!(
            self,
            Self::DynamicAabbTree | Self::SweepAndPrune | Self::Bvh
        )
    }
}

// ---------------------------------------------------------------------------
// SolverType
// ---------------------------------------------------------------------------

/// Constraint / contact solver algorithm.
///
/// All variants are unit-only so this enum can be exposed directly to JS.
#[wasm_bindgen(js_name = "PhysicsSolverType")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SolverType {
    /// Sequential impulses (default for real-time games).
    #[default]
    SequentialImpulse,
    /// Projected Gauss-Seidel (similar to SI, fewer warm-start artefacts).
    ProjectedGaussSeidel,
    /// Position-based dynamics (XPBD). Great for cloth/soft bodies.
    Xpbd,
    /// Non-linear Gauss-Seidel (more accurate for joints).
    NonLinearGaussSeidel,
}

// ---------------------------------------------------------------------------
// GravityPreset
// ---------------------------------------------------------------------------

/// Named gravity presets for common celestial bodies.
///
/// NOTE: This enum has a payload variant `Custom([f64; 3])` and therefore
/// cannot be directly annotated with `#[wasm_bindgen]`. Use
/// `PhysicsConfig::set_gravity_preset_js` or the JSON API to set custom gravity.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum GravityPreset {
    /// Earth surface gravity −9.81 m/s² (default).
    #[default]
    Earth,
    /// Moon surface gravity −1.62 m/s².
    Moon,
    /// Mars surface gravity −3.72 m/s².
    Mars,
    /// Jupiter surface gravity −24.79 m/s².
    Jupiter,
    /// Zero gravity (space / top-down 2D).
    Zero,
    /// Custom gravity vector `[gx, gy, gz]`.
    Custom([f64; 3]),
}

impl GravityPreset {
    /// Return the gravity vector `[gx, gy, gz]` for this preset.
    pub fn to_vec3(&self) -> [f64; 3] {
        match self {
            Self::Earth => [0.0, -9.81, 0.0],
            Self::Moon => [0.0, -1.62, 0.0],
            Self::Mars => [0.0, -3.72, 0.0],
            Self::Jupiter => [0.0, -24.79, 0.0],
            Self::Zero => [0.0, 0.0, 0.0],
            Self::Custom(v) => *v,
        }
    }

    /// Return the magnitude (scalar) of the gravity vector.
    pub fn magnitude(&self) -> f64 {
        let g = self.to_vec3();
        (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
    }
}

// ---------------------------------------------------------------------------
// PhysicsConfig
// ---------------------------------------------------------------------------

/// Comprehensive physics configuration for a browser WASM simulation.
///
/// This type combines all settings from [`SimulationConfig`] with new
/// broadphase, solver, and debugging options.
///
/// Use [`PhysicsConfig::builder()`] to construct with method chaining.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::physics_config::{PhysicsConfig, GravityPreset};
///
/// let cfg = PhysicsConfig::builder()
///     .gravity_preset(GravityPreset::Moon)
///     .substeps(8)
///     .build();
///
/// assert!((cfg.gravity[1] + 1.62).abs() < 1e-6);
/// ```
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicsConfig {
    /// Gravity vector `[gx, gy, gz]` in m/s².
    // `[f64; 3]` does not implement IntoWasmAbi — expose via JS getters.
    #[wasm_bindgen(skip)]
    pub gravity: [f64; 3],
    /// Number of integration substeps per `step()` call.
    pub substeps: u32,
    /// Fixed simulation time step per substep (seconds).
    pub fixed_dt: f64,
    /// Broadphase algorithm.
    pub broadphase: BroadphaseType,
    /// Solver algorithm.
    pub solver: SolverType,
    /// Number of constraint solver iterations per substep.
    pub solver_iterations: u32,
    /// Number of position correction iterations.
    pub position_iterations: u32,
    /// Enable body sleeping when at rest.
    pub sleeping_enabled: bool,
    /// Linear velocity sleep threshold (m/s).
    pub linear_sleep_threshold: f64,
    /// Angular velocity sleep threshold (rad/s).
    pub angular_sleep_threshold: f64,
    /// Time bodies must be slow before sleeping (s).
    pub time_before_sleep: f64,
    /// Enable continuous collision detection (prevents tunneling).
    pub ccd_enabled: bool,
    /// Broadphase AABB margin (m).
    pub broadphase_margin: f64,
    /// Contact penetration slop (m).
    pub contact_slop: f64,
    /// Maximum penetration correction per step (m).
    pub max_penetration_correction: f64,
    /// Enable warm starting of constraints (reuse previous frame impulses).
    pub warm_starting: bool,
    /// Whether debug rendering info should be collected.
    pub debug_rendering: bool,
    /// Target frames-per-second (used by the JS step loop, not physics).
    pub target_fps: u32,
}

impl Default for PhysicsConfig {
    fn default() -> Self {
        let sim = SimulationConfig::default();
        Self {
            gravity: sim.gravity,
            substeps: sim.max_substeps,
            fixed_dt: sim.fixed_dt,
            broadphase: BroadphaseType::default(),
            solver: SolverType::default(),
            solver_iterations: sim.solver_iterations,
            position_iterations: sim.position_iterations,
            sleeping_enabled: sim.sleeping_enabled,
            linear_sleep_threshold: sim.linear_sleep_threshold,
            angular_sleep_threshold: sim.angular_sleep_threshold,
            time_before_sleep: sim.time_before_sleep,
            ccd_enabled: sim.ccd_enabled,
            broadphase_margin: sim.broadphase_margin,
            contact_slop: sim.contact_slop,
            max_penetration_correction: sim.max_penetration_correction,
            warm_starting: sim.warm_starting,
            debug_rendering: false,
            target_fps: 60,
        }
    }
}

impl PhysicsConfig {
    /// Create a default `PhysicsConfig`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Start building a `PhysicsConfig`.
    pub fn builder() -> PhysicsConfigBuilder {
        PhysicsConfigBuilder::new()
    }

    /// Convert to a [`SimulationConfig`] for use with `WasmPhysicsEngine::from_config`.
    pub fn to_simulation_config(&self) -> SimulationConfig {
        SimulationConfig {
            gravity: self.gravity,
            fixed_dt: self.fixed_dt,
            max_substeps: self.substeps,
            solver_iterations: self.solver_iterations,
            linear_sleep_threshold: self.linear_sleep_threshold,
            angular_sleep_threshold: self.angular_sleep_threshold,
            time_before_sleep: self.time_before_sleep,
            sleeping_enabled: self.sleeping_enabled,
            ccd_enabled: self.ccd_enabled,
            broadphase_margin: self.broadphase_margin,
            position_iterations: self.position_iterations,
            warm_starting: self.warm_starting,
            contact_slop: self.contact_slop,
            max_penetration_correction: self.max_penetration_correction,
        }
    }

    /// Serialize to JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }

    /// Deserialize from JSON.
    pub fn from_json(json: &str) -> Option<Self> {
        serde_json::from_str(json).ok()
    }

    /// Return the gravity magnitude (scalar).
    pub fn gravity_magnitude(&self) -> f64 {
        let g = self.gravity;
        (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
    }

    /// Validate this configuration and return a list of warnings.
    ///
    /// Warnings are non-fatal; the config is still usable.
    pub fn validate(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        if self.fixed_dt <= 0.0 {
            warnings.push(format!("fixed_dt must be positive, got {}", self.fixed_dt));
        }
        if self.fixed_dt > 0.1 {
            warnings.push(format!(
                "fixed_dt={} is very large (>0.1 s); simulation may be unstable",
                self.fixed_dt
            ));
        }
        if self.substeps == 0 {
            warnings.push("substeps must be >= 1".to_string());
        }
        if self.solver_iterations == 0 {
            warnings.push("solver_iterations must be >= 1".to_string());
        }
        if self.linear_sleep_threshold < 0.0 {
            warnings.push("linear_sleep_threshold must be non-negative".to_string());
        }
        if self.target_fps == 0 {
            warnings.push("target_fps must be >= 1".to_string());
        }
        warnings
    }
}

// ---------------------------------------------------------------------------
// PhysicsConfig wasm_bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl PhysicsConfig {
    /// Create a default `PhysicsConfig` (accessible from JavaScript).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> PhysicsConfig {
        PhysicsConfig::default()
    }

    /// Gravity vector as `[gx, gy, gz]`.
    #[wasm_bindgen(js_name = "get_gravity")]
    pub fn get_gravity_js(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }

    /// Set gravity directly from three components.
    #[wasm_bindgen(js_name = "set_gravity")]
    pub fn set_gravity_js(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
    }

    /// Set gravity from a named preset string: `"Earth"`, `"Moon"`, `"Mars"`,
    /// `"Jupiter"`, or `"Zero"`. For a custom vector, call `set_gravity` directly.
    #[wasm_bindgen(js_name = "set_gravity_preset")]
    pub fn set_gravity_preset_js(&mut self, preset_name: String) -> Result<(), JsValue> {
        let preset: GravityPreset =
            serde_json::from_str(&format!("\"{}\"", preset_name)).map_err(err_to_jsvalue)?;
        self.gravity = preset.to_vec3();
        Ok(())
    }

    /// Serialize this config to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }

    /// Deserialize a `PhysicsConfig` from a JSON string.
    ///
    /// Returns `null` on parse failure.
    #[wasm_bindgen(js_name = "from_json")]
    pub fn from_json_js(json: String) -> Option<PhysicsConfig> {
        PhysicsConfig::from_json(&json)
    }

    /// Return the magnitude (scalar) of the gravity vector.
    #[wasm_bindgen(js_name = "gravity_magnitude")]
    pub fn gravity_magnitude_js(&self) -> f64 {
        self.gravity_magnitude()
    }

    /// Validate the config and return a `JsValue` array of warning strings.
    ///
    /// Returns an empty array when the config is valid.
    #[wasm_bindgen(js_name = "validate")]
    pub fn validate_js(&self) -> Result<JsValue, JsValue> {
        let warnings = self.validate();
        to_js_value(&warnings)
    }

    /// Convert to a `SimulationConfig` and return as a `JsValue` object.
    #[wasm_bindgen(js_name = "to_simulation_config")]
    pub fn to_simulation_config_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.to_simulation_config())
    }
}

// ---------------------------------------------------------------------------
// PhysicsConfigBuilder
// ---------------------------------------------------------------------------

/// Builder for [`PhysicsConfig`] using method chaining.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::physics_config::{PhysicsConfig, BroadphaseType, GravityPreset};
///
/// let cfg = PhysicsConfig::builder()
///     .gravity_preset(GravityPreset::Moon)
///     .broadphase(BroadphaseType::Bvh)
///     .substeps(8)
///     .ccd(true)
///     .debug_rendering(true)
///     .build();
///
/// assert!(cfg.ccd_enabled);
/// assert!(cfg.debug_rendering);
/// ```
#[derive(Debug, Clone)]
pub struct PhysicsConfigBuilder {
    config: PhysicsConfig,
}

impl PhysicsConfigBuilder {
    /// Create a builder with default values.
    pub fn new() -> Self {
        Self {
            config: PhysicsConfig::default(),
        }
    }

    /// Set gravity from a preset.
    pub fn gravity_preset(mut self, preset: GravityPreset) -> Self {
        self.config.gravity = preset.to_vec3();
        self
    }

    /// Set gravity directly.
    pub fn gravity(mut self, gx: f64, gy: f64, gz: f64) -> Self {
        self.config.gravity = [gx, gy, gz];
        self
    }

    /// Set the number of substeps.
    pub fn substeps(mut self, n: u32) -> Self {
        self.config.substeps = n.max(1);
        self
    }

    /// Set the fixed time step.
    pub fn fixed_dt(mut self, dt: f64) -> Self {
        self.config.fixed_dt = dt;
        self
    }

    /// Set the broadphase algorithm.
    pub fn broadphase(mut self, bp: BroadphaseType) -> Self {
        self.config.broadphase = bp;
        self
    }

    /// Set the solver algorithm.
    pub fn solver(mut self, solver: SolverType) -> Self {
        self.config.solver = solver;
        self
    }

    /// Set the number of solver iterations.
    pub fn solver_iterations(mut self, n: u32) -> Self {
        self.config.solver_iterations = n.max(1);
        self
    }

    /// Enable or disable body sleeping.
    pub fn sleeping(mut self, enabled: bool) -> Self {
        self.config.sleeping_enabled = enabled;
        self
    }

    /// Enable or disable CCD.
    pub fn ccd(mut self, enabled: bool) -> Self {
        self.config.ccd_enabled = enabled;
        self
    }

    /// Enable or disable debug rendering collection.
    pub fn debug_rendering(mut self, enabled: bool) -> Self {
        self.config.debug_rendering = enabled;
        self
    }

    /// Set the target FPS (for the JS step loop).
    pub fn target_fps(mut self, fps: u32) -> Self {
        self.config.target_fps = fps.max(1);
        self
    }

    /// Set the broadphase AABB margin.
    pub fn broadphase_margin(mut self, margin: f64) -> Self {
        self.config.broadphase_margin = margin.max(0.0);
        self
    }

    /// Enable or disable warm starting.
    pub fn warm_starting(mut self, enabled: bool) -> Self {
        self.config.warm_starting = enabled;
        self
    }

    /// Set the contact slop.
    pub fn contact_slop(mut self, slop: f64) -> Self {
        self.config.contact_slop = slop.max(0.0);
        self
    }

    /// Build the `PhysicsConfig`.
    pub fn build(self) -> PhysicsConfig {
        self.config
    }
}

impl Default for PhysicsConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Preset configurations
// ---------------------------------------------------------------------------

/// Return a config tuned for real-time browser physics (mobile-friendly).
pub fn preset_realtime_browser() -> PhysicsConfig {
    PhysicsConfig::builder()
        .gravity_preset(GravityPreset::Earth)
        .substeps(2)
        .fixed_dt(1.0 / 60.0)
        .broadphase(BroadphaseType::SweepAndPrune)
        .solver_iterations(4)
        .sleeping(true)
        .warm_starting(true)
        .build()
}

/// Return a config tuned for high-accuracy offline simulation.
pub fn preset_high_accuracy() -> PhysicsConfig {
    PhysicsConfig::builder()
        .gravity_preset(GravityPreset::Earth)
        .substeps(8)
        .fixed_dt(1.0 / 120.0)
        .broadphase(BroadphaseType::DynamicAabbTree)
        .solver(SolverType::ProjectedGaussSeidel)
        .solver_iterations(16)
        .ccd(true)
        .sleeping(true)
        .build()
}

/// Return a config for a 2D top-down game (no gravity).
pub fn preset_topdown_2d() -> PhysicsConfig {
    PhysicsConfig::builder()
        .gravity_preset(GravityPreset::Zero)
        .substeps(2)
        .fixed_dt(1.0 / 60.0)
        .broadphase(BroadphaseType::UniformGrid)
        .sleeping(true)
        .build()
}

// ---------------------------------------------------------------------------
// Free wasm_bindgen functions for preset access
// ---------------------------------------------------------------------------

/// Create a realtime-browser `PhysicsConfig` from JavaScript.
#[wasm_bindgen(js_name = "physics_config_realtime_browser")]
pub fn physics_config_realtime_browser_js() -> PhysicsConfig {
    preset_realtime_browser()
}

/// Create a high-accuracy `PhysicsConfig` from JavaScript.
#[wasm_bindgen(js_name = "physics_config_high_accuracy")]
pub fn physics_config_high_accuracy_js() -> PhysicsConfig {
    preset_high_accuracy()
}

/// Create a 2D top-down `PhysicsConfig` from JavaScript.
#[wasm_bindgen(js_name = "physics_config_topdown_2d")]
pub fn physics_config_topdown_2d_js() -> PhysicsConfig {
    preset_topdown_2d()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_broadphase_type_name() {
        assert_eq!(BroadphaseType::BruteForce.name(), "BruteForce");
        assert_eq!(BroadphaseType::SweepAndPrune.name(), "SweepAndPrune");
        assert_eq!(BroadphaseType::DynamicAabbTree.name(), "DynamicAabbTree");
    }

    #[test]
    fn test_broadphase_type_dynamic_support() {
        assert!(BroadphaseType::DynamicAabbTree.supports_dynamic_updates());
        assert!(!BroadphaseType::BruteForce.supports_dynamic_updates());
    }

    #[test]
    fn test_gravity_preset_earth() {
        let g = GravityPreset::Earth.to_vec3();
        assert!((g[1] + 9.81).abs() < 1e-10);
        assert!(g[0].abs() < 1e-10);
    }

    #[test]
    fn test_gravity_preset_moon() {
        let g = GravityPreset::Moon.to_vec3();
        assert!((g[1] + 1.62).abs() < 1e-10);
    }

    #[test]
    fn test_gravity_preset_zero() {
        let g = GravityPreset::Zero.to_vec3();
        assert_eq!(g, [0.0, 0.0, 0.0]);
        assert!((GravityPreset::Zero.magnitude()).abs() < 1e-10);
    }

    #[test]
    fn test_gravity_preset_custom() {
        let g = GravityPreset::Custom([1.0, -5.0, 0.0]).to_vec3();
        assert!((g[1] + 5.0).abs() < 1e-10);
        let mag = GravityPreset::Custom([3.0, 4.0, 0.0]).magnitude();
        assert!((mag - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_physics_config_default() {
        let cfg = PhysicsConfig::default();
        assert!((cfg.gravity[1] + 9.81).abs() < 1e-10);
        assert_eq!(cfg.substeps, 4);
        assert_eq!(cfg.target_fps, 60);
    }

    #[test]
    fn test_physics_config_builder_gravity_preset() {
        let cfg = PhysicsConfig::builder()
            .gravity_preset(GravityPreset::Moon)
            .build();
        assert!((cfg.gravity[1] + 1.62).abs() < 1e-10);
    }

    #[test]
    fn test_physics_config_builder_substeps() {
        let cfg = PhysicsConfig::builder().substeps(8).build();
        assert_eq!(cfg.substeps, 8);
    }

    #[test]
    fn test_physics_config_builder_substeps_clamped() {
        let cfg = PhysicsConfig::builder().substeps(0).build();
        assert_eq!(cfg.substeps, 1);
    }

    #[test]
    fn test_physics_config_builder_broadphase() {
        let cfg = PhysicsConfig::builder()
            .broadphase(BroadphaseType::DynamicAabbTree)
            .build();
        assert_eq!(cfg.broadphase, BroadphaseType::DynamicAabbTree);
    }

    #[test]
    fn test_physics_config_builder_ccd() {
        let cfg = PhysicsConfig::builder().ccd(true).build();
        assert!(cfg.ccd_enabled);
    }

    #[test]
    fn test_physics_config_builder_debug_rendering() {
        let cfg = PhysicsConfig::builder().debug_rendering(true).build();
        assert!(cfg.debug_rendering);
    }

    #[test]
    fn test_physics_config_to_simulation_config() {
        let cfg = PhysicsConfig::builder()
            .gravity_preset(GravityPreset::Mars)
            .substeps(6)
            .build();
        let sim = cfg.to_simulation_config();
        assert!((sim.gravity[1] + 3.72).abs() < 1e-10);
        assert_eq!(sim.max_substeps, 6);
    }

    #[test]
    fn test_physics_config_json_roundtrip() {
        let cfg = PhysicsConfig::builder()
            .gravity_preset(GravityPreset::Jupiter)
            .broadphase(BroadphaseType::Bvh)
            .build();
        let json = cfg.to_json();
        let back = PhysicsConfig::from_json(&json).expect("should deserialize");
        assert!((back.gravity[1] + 24.79).abs() < 1e-6);
        assert_eq!(back.broadphase, BroadphaseType::Bvh);
    }

    #[test]
    fn test_physics_config_validate_ok() {
        let cfg = PhysicsConfig::default();
        let warnings = cfg.validate();
        assert!(
            warnings.is_empty(),
            "default config should have no warnings: {:?}",
            warnings
        );
    }

    #[test]
    fn test_physics_config_validate_bad_dt() {
        let cfg = PhysicsConfig {
            fixed_dt: -0.01,
            ..Default::default()
        };
        let warnings = cfg.validate();
        assert!(!warnings.is_empty());
    }

    #[test]
    fn test_physics_config_validate_zero_substeps() {
        let cfg = PhysicsConfig {
            substeps: 0,
            ..Default::default()
        };
        let warnings = cfg.validate();
        assert!(!warnings.is_empty());
    }

    #[test]
    fn test_physics_config_gravity_magnitude() {
        let cfg = PhysicsConfig::builder()
            .gravity_preset(GravityPreset::Earth)
            .build();
        assert!((cfg.gravity_magnitude() - 9.81).abs() < 1e-6);
    }

    #[test]
    fn test_preset_realtime_browser() {
        let cfg = preset_realtime_browser();
        assert_eq!(cfg.substeps, 2);
        assert_eq!(cfg.broadphase, BroadphaseType::SweepAndPrune);
    }

    #[test]
    fn test_preset_high_accuracy() {
        let cfg = preset_high_accuracy();
        assert!(cfg.solver_iterations >= 16);
        assert!(cfg.ccd_enabled);
    }

    #[test]
    fn test_preset_topdown_2d() {
        let cfg = preset_topdown_2d();
        let mag = (cfg.gravity[0] * cfg.gravity[0]
            + cfg.gravity[1] * cfg.gravity[1]
            + cfg.gravity[2] * cfg.gravity[2])
            .sqrt();
        assert!(mag < 1e-10, "top-down 2D preset should have zero gravity");
    }

    #[test]
    fn test_solver_type_default() {
        let cfg = PhysicsConfig::default();
        assert_eq!(cfg.solver, SolverType::SequentialImpulse);
    }

    #[test]
    fn test_physics_config_builder_solver() {
        let cfg = PhysicsConfig::builder().solver(SolverType::Xpbd).build();
        assert_eq!(cfg.solver, SolverType::Xpbd);
    }

    #[test]
    fn test_physics_config_builder_target_fps_clamped() {
        let cfg = PhysicsConfig::builder().target_fps(0).build();
        assert_eq!(cfg.target_fps, 1);
    }

    #[test]
    fn test_physics_config_builder_gravity_direct() {
        let cfg = PhysicsConfig::builder().gravity(0.0, -5.0, 0.0).build();
        assert!((cfg.gravity[1] + 5.0).abs() < 1e-10);
    }

    // ---------------------------------------------------------------------------
    // Integration tests (Slice W6)
    // ---------------------------------------------------------------------------

    /// Verify that the wasm_new constructor produces the same result as default().
    #[test]
    fn test_physics_config_wasm_new_equals_default() {
        let a = PhysicsConfig::default();
        let b = PhysicsConfig::wasm_new();
        assert_eq!(a.substeps, b.substeps);
        assert_eq!(a.target_fps, b.target_fps);
        assert_eq!(a.solver, b.solver);
        assert_eq!(a.broadphase, b.broadphase);
        for i in 0..3 {
            assert!((a.gravity[i] - b.gravity[i]).abs() < 1e-10);
        }
    }

    /// Verify the get_gravity_js / set_gravity_js pair is consistent.
    #[test]
    fn test_get_set_gravity_js_roundtrip() {
        let mut cfg = PhysicsConfig::default();
        cfg.set_gravity_js(1.0, -2.0, 3.0);
        let g = cfg.get_gravity_js();
        assert_eq!(g.len(), 3);
        assert!((g[0] - 1.0).abs() < 1e-10);
        assert!((g[1] + 2.0).abs() < 1e-10);
        assert!((g[2] - 3.0).abs() < 1e-10);
    }
}
