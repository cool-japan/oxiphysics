// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial force-field system for applying continuous forces to rigid bodies.
//!
//! Force fields generate a force vector as a function of world-space position
//! and simulation time.  Multiple fields can be stacked; the
//! `ForceFieldSystem` accumulates them and applies the result to any body
//! that lies within each field's optional bounding region.
//!
//! ## Available field kinds
//!
//! | Kind              | Description                                        |
//! |-------------------|----------------------------------------------------|
//! | `Uniform`       | Constant force independent of position             |
//! | `RadialAttract` | Inverse-power attraction toward a point            |
//! | `RadialRepel`   | Inverse-power repulsion away from a point          |
//! | `Vortex`        | Tangential rotation around an axis                 |
//! | `Wind`          | Directional wind with optional drag scaling        |
//! | `Explosion`     | Expanding pressure wave from an origin point       |
//! | `Turbulent`     | Deterministic pseudo-random perturbation           |
//!
//!
//!
//!
//!
//!
//!
//!
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::force_field::{ForceFieldKind, ForceFieldSystem};
//!
//! let mut sys = ForceFieldSystem::new();
//!
//! // Gravity-like field pulling toward origin
//! sys.add(ForceFieldKind::RadialAttract {
//!     center: [0.0, 0.0, 0.0],
//!     strength: 9.81,
//!     falloff_exp: 2.0,
//! });
//!
//! // Wind pushing in +X
//! sys.add(ForceFieldKind::Wind {
//!     direction: [1.0, 0.0, 0.0],
//!     base_speed: 5.0,
//! });
//!
//! let f = sys.force_at([1.0, 1.0, 1.0], 0.0);
//! assert!(f[0].is_finite()); // wind contribution
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

// ============================================================================
// ForceFieldKind
// ============================================================================

/// Describes the mathematical model of a force field.
#[derive(Debug, Clone)]
pub enum ForceFieldKind {
    /// Constant force applied uniformly at all positions.
    ///
    /// Useful for overriding gravity or applying a global push.
    Uniform {
        /// Force vector in world space \[N\].
        force: [f64; 3],
    },

    /// Gravitational-style attraction toward a point.
    ///
    /// Force magnitude = `strength / r^falloff_exp`.  Clamped so it does not
    /// blow up when `r → 0`.
    RadialAttract {
        /// Attractor position in world space.
        center: [f64; 3],
        /// Strength coefficient (positive → attracts).
        strength: f64,
        /// Exponent for distance falloff (2.0 = inverse-square).
        falloff_exp: f64,
    },

    /// Repulsive force pushing bodies away from a point.
    ///
    /// Force magnitude = `strength / r^falloff_exp`.
    RadialRepel {
        /// Repulsor position in world space.
        center: [f64; 3],
        /// Strength coefficient.
        strength: f64,
        /// Exponent for distance falloff.
        falloff_exp: f64,
    },

    /// Swirling vortex around a line (center + axis).
    ///
    /// Applies a tangential force perpendicular to both the axis and the
    /// radial direction.  Optionally also applies an axial ("suction") force.
    Vortex {
        /// A point on the vortex axis.
        center: [f64; 3],
        /// Unit direction of the vortex axis.
        axis: [f64; 3],
        /// Tangential force per unit radial distance \[N/m\].
        tangential_strength: f64,
        /// Force pulling/pushing along the axis direction \[N\].
        axial_strength: f64,
    },

    /// Directional wind force.
    ///
    /// Force = `mass × drag_coeff × base_speed² × direction` (simplified).
    /// For simplicity this implementation applies force proportional to
    /// `base_speed` and `direction` regardless of the body's own velocity.
    Wind {
        /// Normalised wind direction.
        direction: [f64; 3],
        /// Wind speed \[m/s\] used to compute aerodynamic pressure.
        base_speed: f64,
    },

    /// Spherical pressure wave expanding from a point.
    ///
    /// At simulation time `t`, the wave front is at radius
    /// `wave_speed × (t − start_time)` from the center.  Bodies within
    /// ±`thickness` of the front receive an outward impulse that decays as a
    /// Gaussian.  The field produces zero force before `start_time` and after
    /// the wave has passed.
    Explosion {
        /// Explosion origin in world space.
        center: [f64; 3],
        /// Peak outward force at the wave front \[N\].
        peak_force: f64,
        /// Wave propagation speed \[m/s\].
        wave_speed: f64,
        /// Gaussian half-width of the pressure shell \[m\].
        thickness: f64,
        /// Simulation time when the explosion begins.
        start_time: f64,
    },

    /// Deterministic turbulence — a spatially varying oscillating field.
    ///
    /// Produces a pseudo-random force based on position hashing and a
    /// sinusoidal time modulation.  Useful for wind gusts and jitter.
    Turbulent {
        /// Mean force direction and magnitude.
        base_force: [f64; 3],
        /// Maximum random amplitude added on top of `base_force`.
        amplitude: f64,
        /// Oscillation frequency \[Hz\].
        frequency: f64,
    },
}

impl ForceFieldKind {
    /// Compute the force this field exerts at `pos` at simulation time `t`.
    ///
    /// The result is in Newtons.  Callers should multiply by `1/mass` to get
    /// an acceleration, or multiply by `dt` to get an impulse.
    pub fn force_at(&self, pos: [f64; 3], time: f64) -> [f64; 3] {
        match self {
            ForceFieldKind::Uniform { force } => *force,

            ForceFieldKind::RadialAttract {
                center,
                strength,
                falloff_exp,
            } => {
                let d = vec3_sub(pos, *center);
                let r2 = vec3_dot(d, d);
                if r2 < 1e-12 {
                    return [0.0; 3];
                }
                let r = r2.sqrt();
                let mag = -strength / r.powf(*falloff_exp); // negative = toward center
                let dir = vec3_scale(d, 1.0 / r);
                vec3_scale(dir, mag)
            }

            ForceFieldKind::RadialRepel {
                center,
                strength,
                falloff_exp,
            } => {
                let d = vec3_sub(pos, *center);
                let r2 = vec3_dot(d, d);
                if r2 < 1e-12 {
                    return [0.0; 3];
                }
                let r = r2.sqrt();
                let mag = strength / r.powf(*falloff_exp); // positive = away from center
                let dir = vec3_scale(d, 1.0 / r);
                vec3_scale(dir, mag)
            }

            ForceFieldKind::Vortex {
                center,
                axis,
                tangential_strength,
                axial_strength,
            } => {
                let r = vec3_sub(pos, *center);
                let axis_n = vec3_normalize(*axis);
                let proj_len = vec3_dot(r, axis_n);
                let radial = vec3_sub(r, vec3_scale(axis_n, proj_len));
                let radial_len = vec3_len(radial);
                if radial_len < 1e-12 {
                    return [0.0; 3];
                }
                let tangential_dir = vec3_normalize(vec3_cross(axis_n, radial));
                let tangential = vec3_scale(tangential_dir, tangential_strength * radial_len);
                let axial = vec3_scale(axis_n, *axial_strength);
                vec3_add(tangential, axial)
            }

            ForceFieldKind::Wind {
                direction,
                base_speed,
            } => {
                let dir_n = vec3_normalize(*direction);
                // Simple aerodynamic pressure: F = 0.5 * rho_air * v^2 * Cd * A
                // Here we approximate with F = base_speed * direction (per unit mass)
                vec3_scale(dir_n, *base_speed)
            }

            ForceFieldKind::Explosion {
                center,
                peak_force,
                wave_speed,
                thickness,
                start_time,
            } => {
                if time < *start_time {
                    return [0.0; 3];
                }
                let wave_radius = wave_speed * (time - start_time);
                let d = vec3_sub(pos, *center);
                let r = vec3_len(d);
                if r < 1e-12 {
                    return [0.0; 3];
                }
                let delta = r - wave_radius;
                // Gaussian envelope
                let envelope = (-delta * delta / (2.0 * thickness * thickness)).exp();
                if envelope < 1e-9 {
                    return [0.0; 3];
                }
                let dir = vec3_scale(d, 1.0 / r);
                vec3_scale(dir, peak_force * envelope)
            }

            ForceFieldKind::Turbulent {
                base_force,
                amplitude,
                frequency,
            } => {
                // Deterministic turbulence: hash position into a phase offset,
                // then oscillate at `frequency` Hz.
                let phase_offset = spatial_hash(pos);
                let t_mod = (2.0 * std::f64::consts::PI * frequency * time + phase_offset).sin();
                // Build a perpendicular perturbation to base_force using the hash
                let perturb = [
                    amplitude * t_mod * (phase_offset * 0.123).sin(),
                    amplitude * t_mod * (phase_offset * 0.456).cos(),
                    amplitude * t_mod * (phase_offset * 0.789).sin(),
                ];
                vec3_add(*base_force, perturb)
            }
        }
    }
}

// ============================================================================
// AabbRegion — optional bounding region for a field
// ============================================================================

/// Axis-aligned bounding box that restricts a field to a region.
#[derive(Debug, Clone, Copy)]
pub struct AabbRegion {
    /// Minimum corner \[x, y, z\].
    pub min: [f64; 3],
    /// Maximum corner \[x, y, z\].
    pub max: [f64; 3],
}

impl AabbRegion {
    /// Create a region from min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create a region as a sphere's AABB (conservative test).
    pub fn sphere(center: [f64; 3], radius: f64) -> Self {
        Self {
            min: [center[0] - radius, center[1] - radius, center[2] - radius],
            max: [center[0] + radius, center[1] + radius, center[2] + radius],
        }
    }

    /// Test if `pos` lies strictly inside this region.
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        pos[0] >= self.min[0]
            && pos[0] <= self.max[0]
            && pos[1] >= self.min[1]
            && pos[1] <= self.max[1]
            && pos[2] >= self.min[2]
            && pos[2] <= self.max[2]
    }
}

// ============================================================================
// ForceFieldEntry — one registered field
// ============================================================================

/// A single registered force field with its metadata.
#[derive(Debug, Clone)]
pub struct ForceFieldEntry {
    /// The field's mathematical model.
    pub kind: ForceFieldKind,
    /// Whether this field is currently active.
    pub enabled: bool,
    /// If set, the field only applies to bodies within this AABB.
    pub region: Option<AabbRegion>,
    /// Global scale applied to the field's output force.
    pub multiplier: f64,
}

impl ForceFieldEntry {
    /// Compute the force this entry exerts at `pos` at time `t`.
    ///
    /// Returns `[0,0,0]` if disabled or if `pos` is outside `region`.
    pub fn force_at(&self, pos: [f64; 3], time: f64) -> [f64; 3] {
        if !self.enabled {
            return [0.0; 3];
        }
        if self.region.as_ref().is_some_and(|reg| !reg.contains(pos)) {
            return [0.0; 3];
        }
        vec3_scale(self.kind.force_at(pos, time), self.multiplier)
    }
}

// ============================================================================
// ForceFieldSystem — manages a collection of fields
// ============================================================================

/// Manages a collection of [`ForceFieldEntry`] instances and computes the
/// combined force at any world-space position.
///
/// # Example
///
/// ```rust
/// use oxiphysics::force_field::{ForceFieldKind, ForceFieldSystem, AabbRegion};
///
/// let mut sys = ForceFieldSystem::new();
///
/// // Gravity (downward uniform field)
/// sys.add(ForceFieldKind::Uniform { force: [0.0, -9.81, 0.0] });
///
/// // An explosion at origin, starting at t=0
/// sys.add(ForceFieldKind::Explosion {
///     center: [0.0, 0.0, 0.0],
///     peak_force: 1000.0,
///     wave_speed: 50.0,
///     thickness: 2.0,
///     start_time: 0.0,
/// });
///
/// // Vortex bounded to a region
/// let id = sys.add_bounded(
///     ForceFieldKind::Vortex {
///         center: [0.0, 0.0, 0.0],
///         axis: [0.0, 1.0, 0.0],
///         tangential_strength: 5.0,
///         axial_strength: 0.0,
///     },
///     AabbRegion::new([-5.0; 3], [5.0; 3]),
/// );
///
/// let force = sys.force_at([1.0, 0.0, 0.0], 0.5);
/// assert!(force.iter().any(|&f| f.abs() > 1e-10));
///
/// sys.disable(id);
/// ```
#[derive(Debug, Clone, Default)]
pub struct ForceFieldSystem {
    fields: Vec<ForceFieldEntry>,
}

impl ForceFieldSystem {
    /// Create an empty force field system.
    pub fn new() -> Self {
        Self { fields: Vec::new() }
    }

    /// Add a global (unbounded) force field and return its ID.
    pub fn add(&mut self, kind: ForceFieldKind) -> usize {
        let id = self.fields.len();
        self.fields.push(ForceFieldEntry {
            kind,
            enabled: true,
            region: None,
            multiplier: 1.0,
        });
        id
    }

    /// Add a force field restricted to the given AABB region and return its ID.
    pub fn add_bounded(&mut self, kind: ForceFieldKind, region: AabbRegion) -> usize {
        let id = self.fields.len();
        self.fields.push(ForceFieldEntry {
            kind,
            enabled: true,
            region: Some(region),
            multiplier: 1.0,
        });
        id
    }

    /// Add a field with a custom multiplier.
    pub fn add_scaled(&mut self, kind: ForceFieldKind, multiplier: f64) -> usize {
        let id = self.fields.len();
        self.fields.push(ForceFieldEntry {
            kind,
            enabled: true,
            region: None,
            multiplier,
        });
        id
    }

    /// Disable a field by its ID (preserves the slot; can be re-enabled).
    pub fn disable(&mut self, id: usize) {
        if let Some(f) = self.fields.get_mut(id) {
            f.enabled = false;
        }
    }

    /// Enable a previously disabled field.
    pub fn enable(&mut self, id: usize) {
        if let Some(f) = self.fields.get_mut(id) {
            f.enabled = true;
        }
    }

    /// Remove a field by its ID.  Subsequent IDs are **not** renumbered.
    /// The slot is replaced with a disabled dummy field.
    pub fn remove(&mut self, id: usize) {
        if let Some(f) = self.fields.get_mut(id) {
            f.enabled = false;
            f.multiplier = 0.0;
        }
    }

    /// Set the multiplier for a field (can be used to fade a field in/out).
    pub fn set_multiplier(&mut self, id: usize, multiplier: f64) {
        if let Some(f) = self.fields.get_mut(id) {
            f.multiplier = multiplier;
        }
    }

    /// Compute the total force at `pos` at simulation time `time`.
    ///
    /// Accumulates contributions from all enabled fields (respecting regions).
    pub fn force_at(&self, pos: [f64; 3], time: f64) -> [f64; 3] {
        let mut total = [0.0_f64; 3];
        for f in &self.fields {
            let fv = f.force_at(pos, time);
            total[0] += fv[0];
            total[1] += fv[1];
            total[2] += fv[2];
        }
        total
    }

    /// Apply forces to a batch of (position, force-accumulator) pairs.
    ///
    /// `items` is a slice of `(position, &mut force_accumulator)` where the
    /// force accumulator is modified in-place.
    pub fn apply_to_batch(&self, items: &mut [([f64; 3], [f64; 3])], time: f64) {
        for (pos, acc) in items.iter_mut() {
            let fv = self.force_at(*pos, time);
            acc[0] += fv[0];
            acc[1] += fv[1];
            acc[2] += fv[2];
        }
    }

    /// Number of registered fields (including disabled/removed ones).
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// Returns `true` if no fields are registered.
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// Number of currently enabled fields.
    pub fn active_count(&self) -> usize {
        self.fields
            .iter()
            .filter(|f| f.enabled && f.multiplier.abs() > 1e-15)
            .count()
    }

    /// Read-only access to all registered fields.
    pub fn fields(&self) -> &[ForceFieldEntry] {
        &self.fields
    }
}

// ============================================================================
// Internal vector math
// ============================================================================

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn vec3_len(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a);
    if l > 1e-15 {
        vec3_scale(a, 1.0 / l)
    } else {
        [0.0, 1.0, 0.0]
    }
}

/// Deterministic spatial hash → phase angle in [0, 2π).
#[inline]
fn spatial_hash(pos: [f64; 3]) -> f64 {
    // Use integer hashing on quantized coordinates
    let ix = (pos[0] * 1000.0) as i64;
    let iy = (pos[1] * 1000.0) as i64;
    let iz = (pos[2] * 1000.0) as i64;
    let h = ix.wrapping_mul(2654435761_i64)
        ^ iy.wrapping_mul(1234567891_i64)
        ^ iz.wrapping_mul(9876543211_i64);
    let h_f = (h.unsigned_abs() % 1_000_000) as f64 / 1_000_000.0;
    h_f * 2.0 * std::f64::consts::PI
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_field() {
        let k = ForceFieldKind::Uniform {
            force: [0.0, -9.81, 0.0],
        };
        let f = k.force_at([5.0, 10.0, -3.0], 1.0);
        assert!((f[1] + 9.81).abs() < 1e-12);
        assert!(f[0].abs() < 1e-12);
    }

    #[test]
    fn radial_attract_toward_center() {
        let k = ForceFieldKind::RadialAttract {
            center: [0.0; 3],
            strength: 10.0,
            falloff_exp: 2.0,
        };
        let f = k.force_at([1.0, 0.0, 0.0], 0.0);
        // Force should point in -X direction
        assert!(f[0] < 0.0);
    }

    #[test]
    fn radial_repel_away_from_center() {
        let k = ForceFieldKind::RadialRepel {
            center: [0.0; 3],
            strength: 10.0,
            falloff_exp: 2.0,
        };
        let f = k.force_at([1.0, 0.0, 0.0], 0.0);
        assert!(f[0] > 0.0); // away from center
    }

    #[test]
    fn vortex_perpendicular_to_axis_and_radial() {
        let k = ForceFieldKind::Vortex {
            center: [0.0; 3],
            axis: [0.0, 1.0, 0.0],
            tangential_strength: 5.0,
            axial_strength: 0.0,
        };
        let f = k.force_at([1.0, 0.0, 0.0], 0.0);
        // Force should be in the Z direction (tangent to circle in XZ plane around Y)
        assert!(f[2].abs() > 1e-10);
        assert!(f[1].abs() < 1e-10); // no axial component
    }

    #[test]
    fn explosion_before_start_is_zero() {
        let k = ForceFieldKind::Explosion {
            center: [0.0; 3],
            peak_force: 1000.0,
            wave_speed: 10.0,
            thickness: 1.0,
            start_time: 5.0,
        };
        let f = k.force_at([1.0, 0.0, 0.0], 3.0); // before start_time
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn explosion_peak_at_wave_front() {
        let peak = 1000.0_f64;
        let speed = 10.0_f64;
        let t = 1.0_f64;
        // Place test point exactly at wave front (r = speed * t = 10)
        let k = ForceFieldKind::Explosion {
            center: [0.0; 3],
            peak_force: peak,
            wave_speed: speed,
            thickness: 0.5,
            start_time: 0.0,
        };
        let f = k.force_at([speed * t, 0.0, 0.0], t);
        // Magnitude should be close to peak_force (Gaussian envelope ≈ 1 at wave front)
        let mag = vec3_len(f);
        assert!((mag - peak).abs() < 1.0, "mag={mag}, expected≈{peak}");
    }

    #[test]
    fn turbulent_deterministic() {
        let k = ForceFieldKind::Turbulent {
            base_force: [1.0, 0.0, 0.0],
            amplitude: 0.1,
            frequency: 1.0,
        };
        let f1 = k.force_at([1.0, 2.0, 3.0], 0.5);
        let f2 = k.force_at([1.0, 2.0, 3.0], 0.5);
        assert_eq!(f1, f2); // deterministic
    }

    #[test]
    fn system_accumulates_fields() {
        let mut sys = ForceFieldSystem::new();
        sys.add(ForceFieldKind::Uniform {
            force: [1.0, 0.0, 0.0],
        });
        sys.add(ForceFieldKind::Uniform {
            force: [0.0, 1.0, 0.0],
        });
        let f = sys.force_at([0.0; 3], 0.0);
        assert!((f[0] - 1.0).abs() < 1e-12);
        assert!((f[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn system_disable_field() {
        let mut sys = ForceFieldSystem::new();
        let id = sys.add(ForceFieldKind::Uniform {
            force: [100.0, 0.0, 0.0],
        });
        sys.disable(id);
        let f = sys.force_at([0.0; 3], 0.0);
        assert!(f[0].abs() < 1e-12);
    }

    #[test]
    fn bounded_field_outside_region() {
        let mut sys = ForceFieldSystem::new();
        sys.add_bounded(
            ForceFieldKind::Uniform {
                force: [100.0, 0.0, 0.0],
            },
            AabbRegion::new([0.0; 3], [1.0; 3]),
        );
        // Outside the region — force should be zero
        let f = sys.force_at([5.0, 0.0, 0.0], 0.0);
        assert!(f[0].abs() < 1e-12);
    }

    #[test]
    fn bounded_field_inside_region() {
        let mut sys = ForceFieldSystem::new();
        sys.add_bounded(
            ForceFieldKind::Uniform {
                force: [100.0, 0.0, 0.0],
            },
            AabbRegion::new([0.0; 3], [10.0; 3]),
        );
        let f = sys.force_at([5.0, 5.0, 5.0], 0.0);
        assert!((f[0] - 100.0).abs() < 1e-12);
    }

    #[test]
    fn apply_to_batch() {
        let mut sys = ForceFieldSystem::new();
        sys.add(ForceFieldKind::Uniform {
            force: [1.0, 2.0, 3.0],
        });
        let mut batch = [
            ([0.0_f64; 3], [0.0_f64; 3]),
            ([1.0, 1.0, 1.0], [0.0_f64; 3]),
        ];
        sys.apply_to_batch(&mut batch, 0.0);
        for (_, acc) in &batch {
            assert!((acc[0] - 1.0).abs() < 1e-12);
            assert!((acc[1] - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn active_count() {
        let mut sys = ForceFieldSystem::new();
        let id0 = sys.add(ForceFieldKind::Uniform { force: [0.0; 3] });
        let _id1 = sys.add(ForceFieldKind::Uniform { force: [0.0; 3] });
        assert_eq!(sys.active_count(), 2);
        sys.disable(id0);
        assert_eq!(sys.active_count(), 1);
    }
}
