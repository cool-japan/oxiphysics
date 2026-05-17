// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WASM boundary types: `WasmVec3`, `WasmTransform`, error helpers, memory
//! helpers, and TypeScript type definition constants.

use wasm_bindgen::prelude::*;

// ===========================================================================
// WasmVec3 — JS-facing 3D vector with named getters
// ===========================================================================

/// A lightweight, JS-friendly 3D vector type.
///
/// Mirrors the structure of `Vec3Wasm` but is intended as the canonical
/// "hand-off" type at the wasm-bindgen boundary.
///
/// # Example (native)
///
/// ```no_run
/// use oxiphysics_wasm::engine::WasmVec3;
///
/// let v = WasmVec3::new(1.0, 2.0, 3.0);
/// assert!((v.x() - 1.0).abs() < 1e-12);
/// assert!((v.length() - f64::sqrt(14.0)).abs() < 1e-10);
/// ```
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WasmVec3 {
    x: f64,
    y: f64,
    z: f64,
}

#[wasm_bindgen]
impl WasmVec3 {
    /// Construct from components.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The zero vector.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// X component getter.
    #[wasm_bindgen(getter)]
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Y component getter.
    #[wasm_bindgen(getter)]
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Z component getter.
    #[wasm_bindgen(getter)]
    pub fn z(&self) -> f64 {
        self.z
    }

    /// X component setter.
    #[wasm_bindgen(setter)]
    pub fn set_x(&mut self, v: f64) {
        self.x = v;
    }

    /// Y component setter.
    #[wasm_bindgen(setter)]
    pub fn set_y(&mut self, v: f64) {
        self.y = v;
    }

    /// Z component setter.
    #[wasm_bindgen(setter)]
    pub fn set_z(&mut self, v: f64) {
        self.z = v;
    }

    /// Return as flat `Vec<f64>` of `[x, y, z]` (JS-compatible).
    pub fn to_array_js(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }

    /// Return as `js_sys::Float64Array` (typed array, JS-compatible).
    pub fn to_typed_array(&self) -> js_sys::Float64Array {
        js_sys::Float64Array::from([self.x, self.y, self.z].as_slice())
    }

    /// Euclidean length.
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Return normalised copy, or zero if length is tiny.
    pub fn normalized(&self) -> Self {
        let l = self.length();
        if l < 1e-15 {
            Self::zero()
        } else {
            Self {
                x: self.x / l,
                y: self.y / l,
                z: self.z / l,
            }
        }
    }

    /// Dot product with another `WasmVec3`.
    pub fn dot(&self, other: &WasmVec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Component-wise add.
    pub fn add(&self, other: &WasmVec3) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    /// Component-wise subtract.
    pub fn sub(&self, other: &WasmVec3) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    /// Scalar multiply.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    /// Distance to another vector.
    pub fn distance_to(&self, other: &WasmVec3) -> f64 {
        self.sub(other).length()
    }

    /// Serialize to JSON string (for JS postMessage).
    pub fn to_json(&self) -> String {
        format!(r#"{{"x":{},"y":{},"z":{}}}"#, self.x, self.y, self.z)
    }
}

impl WasmVec3 {
    /// Return as flat `[x, y, z]` array (Rust-only).
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Create from a flat `[f64; 3]` array (Rust-only).
    pub fn from_array(a: [f64; 3]) -> Self {
        Self {
            x: a[0],
            y: a[1],
            z: a[2],
        }
    }
}

impl Default for WasmVec3 {
    fn default() -> Self {
        Self::zero()
    }
}

// ===========================================================================
// WasmTransform — position + orientation as flat JS arrays
// ===========================================================================

/// A rigid-body transform exposed at the JS boundary.
///
/// Position is stored as `(px, py, pz)` and rotation as a `[qx, qy, qz, qw]`
/// unit quaternion. The whole transform can be flattened to a `Vec<f64>` of
/// 7 elements for zero-copy `Float64Array` transfer.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::engine::WasmTransform;
///
/// let t = WasmTransform::identity();
/// let flat = t.to_flat();
/// assert_eq!(flat.len(), 7);
/// assert!((flat[6] - 1.0).abs() < 1e-12); // w == 1 for identity
/// ```
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WasmTransform {
    /// Position X.
    pub px: f64,
    /// Position Y.
    pub py: f64,
    /// Position Z.
    pub pz: f64,
    /// Rotation quaternion X.
    pub qx: f64,
    /// Rotation quaternion Y.
    pub qy: f64,
    /// Rotation quaternion Z.
    pub qz: f64,
    /// Rotation quaternion W.
    pub qw: f64,
}

#[wasm_bindgen]
impl WasmTransform {
    /// Construct from position (px, py, pz) and quaternion (qx, qy, qz, qw).
    pub fn new(px: f64, py: f64, pz: f64, qx: f64, qy: f64, qz: f64, qw: f64) -> Self {
        Self {
            px,
            py,
            pz,
            qx,
            qy,
            qz,
            qw,
        }
    }

    /// Identity transform at origin.
    pub fn identity() -> Self {
        Self {
            px: 0.0,
            py: 0.0,
            pz: 0.0,
            qx: 0.0,
            qy: 0.0,
            qz: 0.0,
            qw: 1.0,
        }
    }

    /// Create from a position with identity rotation.
    pub fn from_position(x: f64, y: f64, z: f64) -> Self {
        Self {
            px: x,
            py: y,
            pz: z,
            qx: 0.0,
            qy: 0.0,
            qz: 0.0,
            qw: 1.0,
        }
    }

    /// Get position as `WasmVec3`.
    pub fn position(&self) -> WasmVec3 {
        WasmVec3::new(self.px, self.py, self.pz)
    }

    /// Get position as `Vec<f64>` of `[x, y, z]` (JS-compatible).
    pub fn position_js(&self) -> Vec<f64> {
        vec![self.px, self.py, self.pz]
    }

    /// Get rotation quaternion as `Vec<f64>` of `[qx, qy, qz, qw]` (JS-compatible).
    pub fn rotation_js(&self) -> Vec<f64> {
        vec![self.qx, self.qy, self.qz, self.qw]
    }

    /// Flatten to `[px, py, pz, qx, qy, qz, qw]` — a `Vec<f64>` suitable
    /// for returning as a JS `Float64Array`.
    pub fn to_flat(&self) -> Vec<f64> {
        vec![
            self.px, self.py, self.pz, self.qx, self.qy, self.qz, self.qw,
        ]
    }

    /// Reconstruct from a flat `[px, py, pz, qx, qy, qz, qw]` `Vec<f64>`.
    ///
    /// Returns identity if the slice has fewer than 7 elements.
    pub fn from_flat_js(data: Vec<f64>) -> Self {
        if data.len() < 7 {
            return Self::identity();
        }
        Self {
            px: data[0],
            py: data[1],
            pz: data[2],
            qx: data[3],
            qy: data[4],
            qz: data[5],
            qw: data[6],
        }
    }

    /// Column-major 4×4 matrix as `Vec<f64>` (16 elements).
    ///
    /// Compatible with WebGL `uniformMatrix4fv`.
    pub fn to_matrix4_js(&self) -> Vec<f64> {
        let (qx, qy, qz, qw) = (self.qx, self.qy, self.qz, self.qw);
        let (tx, ty, tz) = (self.px, self.py, self.pz);
        vec![
            1.0 - 2.0 * (qy * qy + qz * qz),
            2.0 * (qx * qy + qw * qz),
            2.0 * (qx * qz - qw * qy),
            0.0,
            2.0 * (qx * qy - qw * qz),
            1.0 - 2.0 * (qx * qx + qz * qz),
            2.0 * (qy * qz + qw * qx),
            0.0,
            2.0 * (qx * qz + qw * qy),
            2.0 * (qy * qz - qw * qx),
            1.0 - 2.0 * (qx * qx + qy * qy),
            0.0,
            tx,
            ty,
            tz,
            1.0,
        ]
    }

    /// Serialize to JSON string.
    pub fn to_json(&self) -> String {
        format!(
            r#"{{"position":{{"x":{},"y":{},"z":{}}},"rotation":{{"x":{},"y":{},"z":{},"w":{}}}}}"#,
            self.px, self.py, self.pz, self.qx, self.qy, self.qz, self.qw,
        )
    }
}

impl WasmTransform {
    /// Construct from a `WasmVec3` position and `[f64; 4]` quaternion (Rust-only).
    pub fn from_vec3_rot(position: WasmVec3, rotation: [f64; 4]) -> Self {
        Self {
            px: position.x,
            py: position.y,
            pz: position.z,
            qx: rotation[0],
            qy: rotation[1],
            qz: rotation[2],
            qw: rotation[3],
        }
    }

    /// Return position as `WasmVec3` (Rust-only alias).
    pub fn position_vec3(&self) -> WasmVec3 {
        WasmVec3::new(self.px, self.py, self.pz)
    }

    /// Return rotation as `[f64; 4]` (Rust-only).
    pub fn rotation(&self) -> [f64; 4] {
        [self.qx, self.qy, self.qz, self.qw]
    }

    /// Return position as `[f64; 3]` (Rust-only).
    pub fn position_array(&self) -> [f64; 3] {
        [self.px, self.py, self.pz]
    }

    /// Column-major 4×4 matrix as `[f64; 16]` (Rust-only).
    pub fn to_matrix4(&self) -> [f64; 16] {
        let m = self.to_matrix4_js();
        [
            m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13],
            m[14], m[15],
        ]
    }

    /// Reconstruct from a flat `[px, py, pz, qx, qy, qz, qw]` slice (Rust-only).
    ///
    /// Returns `None` if the slice has fewer than 7 elements.
    pub fn from_flat(data: &[f64]) -> Option<Self> {
        if data.len() < 7 {
            return None;
        }
        Some(Self {
            px: data[0],
            py: data[1],
            pz: data[2],
            qx: data[3],
            qy: data[4],
            qz: data[5],
            qw: data[6],
        })
    }
}

impl Default for WasmTransform {
    fn default() -> Self {
        Self::identity()
    }
}

// ===========================================================================
// Error helpers — convert Rust errors to JS-friendly strings / JsValue
// ===========================================================================

/// Convert a physics [`crate::error::Error`] to a JSON-encoded JS error object.
///
/// The resulting string can be passed directly to `new Error(msg)` in JS or
/// thrown as a string from a wasm-bindgen function.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::engine::error_to_js_string;
/// use oxiphysics_wasm::error::Error;
///
/// let err = Error::InvalidHandle(42);
/// let s = error_to_js_string(&format!("{}", err));
/// assert!(s.contains("InvalidHandle"));
/// assert!(s.contains("42"));
/// ```
#[wasm_bindgen]
pub fn error_to_js_string(err_msg: &str) -> String {
    format!(r#"{{"error":"{}"}}"#, err_msg)
}

/// Convert a physics error to a `JsValue` for use in `Result<T, JsValue>`.
pub fn err_to_jsvalue(err: &crate::error::Error) -> JsValue {
    JsValue::from_str(&format!("{err}"))
}

/// Convert a `Result<T, Error>` to `Result<T, JsValue>` for wasm-bindgen
/// functions that return `Result<_, JsValue>`.
pub fn result_to_js<T>(r: crate::error::Result<T>) -> std::result::Result<T, JsValue> {
    r.map_err(|e| err_to_jsvalue(&e))
}

// ===========================================================================
// Memory management helpers — free allocated flat buffers
// ===========================================================================

/// Helper to drop (free) a heap-allocated `Vec<f64>` from WASM.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::engine::free_f64_buffer;
///
/// let buf: Vec<f64> = vec![1.0, 2.0, 3.0];
/// free_f64_buffer(buf); // drops the allocation
/// ```
#[wasm_bindgen]
pub fn free_f64_buffer(buf: Vec<f64>) {
    drop(buf);
}

/// Helper to drop a heap-allocated `Vec<u32>`.
#[wasm_bindgen]
pub fn free_u32_buffer(buf: Vec<u32>) {
    drop(buf);
}

/// Helper to drop a heap-allocated `Vec<u8>` (e.g. serialized JSON).
#[wasm_bindgen]
pub fn free_u8_buffer(buf: Vec<u8>) {
    drop(buf);
}

/// Return the size of the WASM memory page in bytes (64 KiB per WebAssembly spec).
pub const WASM_PAGE_SIZE: usize = 65536;

// ===========================================================================
// TypeScript type definition string constants
// ===========================================================================

/// TypeScript type definitions for the WASM physics engine.
///
/// Embed these in your build toolchain or write them to a `.d.ts` file so
/// TypeScript consumers of the WASM module get full type information.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::engine::TS_DEFINITIONS;
/// assert!(TS_DEFINITIONS.contains("WasmPhysicsEngine"));
/// assert!(TS_DEFINITIONS.contains("bodyCount()"));
/// ```
pub const TS_DEFINITIONS: &str = r#"
// Auto-generated TypeScript definitions for oxiphysics-wasm.
// DO NOT EDIT — regenerate via the Rust build script.

export interface Vec3 {
  x: number;
  y: number;
  z: number;
}

export interface Quaternion {
  x: number;
  y: number;
  z: number;
  w: number;
}

export interface Transform {
  position: Vec3;
  rotation: Quaternion;
}

export interface BodyState {
  handle: number;
  position: Vec3;
  rotation: Quaternion;
  linearVelocity: Vec3;
  angularVelocity: Vec3;
  isSleeping: boolean;
  isActive: boolean;
  kineticEnergy: number;
}

export interface ContactInfo {
  bodyA: number;
  bodyB: number;
  pointOnA: Vec3;
  pointOnB: Vec3;
  normal: Vec3;
  depth: number;
  relativeVelocity: number;
  impulse: number;
  isNew: boolean;
  frictionImpulse: number;
}

export interface RaycastResult {
  hit: boolean;
  bodyHandle: number;
  point: Vec3;
  normal: Vec3;
  distance: number;
}

export interface SimulationConfig {
  gravityX: number;
  gravityY: number;
  gravityZ: number;
  fixedDt: number;
  maxSubsteps: number;
  solverIterations: number;
  ccdEnabled: boolean;
  sleepingEnabled: boolean;
}

/** Main WASM physics engine. */
export class WasmPhysicsEngine {
  /** Create engine with given gravity. */
  static new(gx: number, gy: number, gz: number): WasmPhysicsEngine;

  /** Add a rigid body and return its handle. */
  addRigidBody(x: number, y: number, z: number, mass: number): number;

  /** Set linear velocity of body with given handle. */
  setVelocity(id: number, vx: number, vy: number, vz: number): void;

  /** Get position as [x, y, z]. */
  getPosition(id: number): Float64Array;

  /** Advance simulation by dt seconds. */
  step(dt: number): void;

  /** Get all body positions as flat Float64Array [x0,y0,z0, x1,y1,z1, ...]. */
  getAllPositions(): Float64Array;

  /** Set gravity vector. */
  setGravity(gx: number, gy: number, gz: number): void;

  /** Return number of active bodies. */
  bodyCount(): number;

  /** Add a sphere collider to a body. */
  addColliderSphere(bodyId: number, radius: number): number;

  /** Add a box collider to a body. */
  addColliderBox(bodyId: number, hx: number, hy: number, hz: number): number;

  /** Get contact information from last step. */
  getContacts(): ContactInfo[];

  /** Run a WebGPU compute placeholder (returns mock data). */
  webgpuComputePlaceholder(workgroupSize: number): Float64Array;
}

/** JS-facing 3D vector. */
export class WasmVec3 {
  constructor(x: number, y: number, z: number);
  readonly x: number;
  readonly y: number;
  readonly z: number;
  length(): number;
  normalized(): WasmVec3;
  dot(other: WasmVec3): number;
  toArrayJs(): Float64Array;
  toJson(): string;
}

/** JS-facing rigid-body transform. */
export class WasmTransform {
  static identity(): WasmTransform;
  static fromPosition(x: number, y: number, z: number): WasmTransform;
  position(): WasmVec3;
  positionJs(): Float64Array;
  rotationJs(): Float64Array;
  toFlat(): Float64Array;
  toMatrix4Js(): Float64Array;
  toJson(): string;
}
"#;

/// Minimal TypeScript definitions for the Vec3 type only.
pub const TS_VEC3_DEF: &str = r#"
export interface Vec3 { x: number; y: number; z: number; }
export declare function vec3(x: number, y: number, z: number): Vec3;
"#;

/// Minimal TypeScript definitions for Transform.
pub const TS_TRANSFORM_DEF: &str = r#"
export interface Transform {
  position: { x: number; y: number; z: number };
  rotation: { x: number; y: number; z: number; w: number };
  toFlat(): Float64Array;
  toMatrix4Js(): Float64Array;
}
"#;

/// Minimal TypeScript definitions for contact info.
pub const TS_CONTACT_DEF: &str = r#"
export interface ContactInfo {
  bodyA: number;
  bodyB: number;
  normal: { x: number; y: number; z: number };
  depth: number;
  impulse: number;
}
"#;
