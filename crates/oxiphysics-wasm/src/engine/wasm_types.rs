// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WASM boundary types: `WasmVec3`, `WasmTransform`, error helpers, memory
//! helpers, and TypeScript type definition constants.

#![allow(missing_docs)]

// ===========================================================================
// WasmVec3 — JS-facing 3D vector with named getters
// ===========================================================================

/// A lightweight, JS-friendly 3D vector type.
///
/// Mirrors the structure of `Vec3Wasm` but is intended as the canonical
/// "hand-off" type at the wasm-bindgen boundary.  In a real wasm-bindgen build
/// this struct would carry `#[wasm_bindgen]`; here we keep it as plain Rust so
/// that native tests work without a wasm toolchain.
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
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub struct WasmVec3 {
    x: f64,
    y: f64,
    z: f64,
}

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

    /// X component getter (mirrors wasm-bindgen `#[wasm_bindgen(getter)]`).
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Y component getter.
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Z component getter.
    pub fn z(&self) -> f64 {
        self.z
    }

    /// X component setter.
    pub fn set_x(&mut self, v: f64) {
        self.x = v;
    }

    /// Y component setter.
    pub fn set_y(&mut self, v: f64) {
        self.y = v;
    }

    /// Z component setter.
    pub fn set_z(&mut self, v: f64) {
        self.z = v;
    }

    /// Return as flat `[x, y, z]` array.
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Create from a flat `[f64; 3]` array.
    pub fn from_array(a: [f64; 3]) -> Self {
        Self {
            x: a[0],
            y: a[1],
            z: a[2],
        }
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

    /// Dot product.
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
/// Position is a `WasmVec3` and rotation is a `[x, y, z, w]` unit quaternion.
/// The whole transform can be flattened to a `[f64; 7]` for zero-copy
/// `Float64Array` transfer.
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
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub struct WasmTransform {
    /// Position.
    position: WasmVec3,
    /// Rotation quaternion `[x, y, z, w]`.
    rotation: [f64; 4],
}

impl WasmTransform {
    /// Construct from position and rotation.
    pub fn new(position: WasmVec3, rotation: [f64; 4]) -> Self {
        Self { position, rotation }
    }

    /// Identity transform at origin.
    pub fn identity() -> Self {
        Self {
            position: WasmVec3::zero(),
            rotation: [0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Create from a position with identity rotation.
    pub fn from_position(x: f64, y: f64, z: f64) -> Self {
        Self {
            position: WasmVec3::new(x, y, z),
            rotation: [0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Position getter.
    pub fn position(&self) -> WasmVec3 {
        self.position
    }

    /// Rotation quaternion getter `[qx, qy, qz, qw]`.
    pub fn rotation(&self) -> [f64; 4] {
        self.rotation
    }

    /// Return position as flat `[f64; 3]`.
    pub fn position_array(&self) -> [f64; 3] {
        self.position.to_array()
    }

    /// Flatten to `[px, py, pz, qx, qy, qz, qw]` — a `Vec`f64` suitable
    /// for returning as a JS `Float64Array`.
    pub fn to_flat(&self) -> Vec<f64> {
        let p = self.position.to_array();
        let q = self.rotation;
        vec![p[0], p[1], p[2], q[0], q[1], q[2], q[3]]
    }

    /// Reconstruct from a flat `\[px, py, pz, qx, qy, qz, qw\]` slice.
    ///
    /// Returns `None` if the slice has fewer than 7 elements.
    pub fn from_flat(data: &[f64]) -> Option<Self> {
        if data.len() < 7 {
            return None;
        }
        Some(Self {
            position: WasmVec3::new(data[0], data[1], data[2]),
            rotation: [data[3], data[4], data[5], data[6]],
        })
    }

    /// Column-major 4×4 matrix representation (16 floats).
    ///
    /// Compatible with WebGL `uniformMatrix4fv`.
    pub fn to_matrix4(&self) -> [f64; 16] {
        let [qx, qy, qz, qw] = self.rotation;
        let [tx, ty, tz] = self.position.to_array();
        [
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
        let [qx, qy, qz, qw] = self.rotation;
        let p = self.position;
        format!(
            r#"{{"position":{{"x":{},"y":{},"z":{}}},"rotation":{{"x":{},"y":{},"z":{},"w":{}}}}}"#,
            p.x(),
            p.y(),
            p.z(),
            qx,
            qy,
            qz,
            qw
        )
    }
}

impl Default for WasmTransform {
    fn default() -> Self {
        Self::identity()
    }
}

// ===========================================================================
// Error helpers — convert Rust errors to JS-friendly strings
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
/// let s = error_to_js_string(&Error::InvalidHandle(42));
/// assert!(s.contains("InvalidHandle"));
/// assert!(s.contains("42"));
/// ```
#[allow(dead_code)]
pub fn error_to_js_string(err: &crate::error::Error) -> String {
    format!(r#"{{"error":"{}","detail":"{}"}}"#, err, err.to_json())
}

/// Convert a `Result<T, Error>` to `Result<T, String>` for wasm-bindgen
/// functions that return `Result<_, JsValue>`.  On the native side the
/// `String` plays the role of `JsValue`.
#[allow(dead_code)]
pub fn result_to_js<T>(r: crate::error::Result<T>) -> std::result::Result<T, String> {
    r.map_err(|e| error_to_js_string(&e))
}

// ===========================================================================
// Memory management helpers — free allocated flat buffers
// ===========================================================================

/// Helper to drop (free) a heap-allocated `Vec`f64` from WASM.
///
/// In a true wasm-bindgen build the GC/FinalizationRegistry handles this;
/// in native tests we just drop the vector explicitly.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::engine::free_f64_buffer;
///
/// let buf: Vec<f64> = vec![1.0, 2.0, 3.0];
/// free_f64_buffer(buf); // drops the allocation
/// ```
#[allow(dead_code)]
pub fn free_f64_buffer(buf: Vec<f64>) {
    drop(buf);
}

/// Helper to drop a heap-allocated `Vec`u32`.
#[allow(dead_code)]
pub fn free_u32_buffer(buf: Vec<u32>) {
    drop(buf);
}

/// Helper to drop a heap-allocated `Vec`u8` (e.g. serialized JSON).
#[allow(dead_code)]
pub fn free_u8_buffer(buf: Vec<u8>) {
    drop(buf);
}

/// Return the size of the WASM memory page in bytes (64 KiB per WebAssembly spec).
#[allow(dead_code)]
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
  toArray(): Float64Array;
  toJson(): string;
}

/** JS-facing rigid-body transform. */
export class WasmTransform {
  static identity(): WasmTransform;
  static fromPosition(x: number, y: number, z: number): WasmTransform;
  readonly position: WasmVec3;
  readonly rotation: Float64Array;
  toFlat(): Float64Array;
  toMatrix4(): Float64Array;
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
  toMatrix4(): Float64Array;
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
