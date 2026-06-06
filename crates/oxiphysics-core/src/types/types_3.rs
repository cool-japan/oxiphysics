//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use crate::math::{Mat3, Quat, Real, Vec3};

/// Mass properties of a rigid body.
#[derive(Debug, Clone)]
pub struct MassProperties {
    /// Mass in kilograms.
    pub mass: Real,
    /// Center of mass in local space.
    pub center_of_mass: Vec3,
    /// Local inertia tensor (3x3 matrix).
    pub local_inertia: Mat3,
}
impl MassProperties {
    /// Create new mass properties.
    pub fn new(mass: Real, center_of_mass: Vec3, local_inertia: Mat3) -> Self {
        Self {
            mass,
            center_of_mass,
            local_inertia,
        }
    }
    /// Mass properties for a point mass.
    pub fn point_mass(mass: Real) -> Self {
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::zeros(),
        }
    }
    /// Inverse mass (0 for infinite/static mass).
    pub fn inverse_mass(&self) -> Real {
        if self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        }
    }
    /// Inverse inertia tensor (zero matrix for infinite/static inertia).
    pub fn inverse_inertia(&self) -> Mat3 {
        if self.mass > 0.0 {
            self.local_inertia.try_inverse().unwrap_or_else(Mat3::zeros)
        } else {
            Mat3::zeros()
        }
    }
}
impl MassProperties {
    /// Mass properties for a uniform-density sphere.
    pub fn sphere(mass: Real, radius: Real) -> Self {
        let i_scalar = (2.0 / 5.0) * mass * radius * radius;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(i_scalar, 0.0, 0.0, 0.0, i_scalar, 0.0, 0.0, 0.0, i_scalar),
        }
    }
    /// Mass properties for a uniform-density box with dimensions (wx, wy, wz).
    pub fn cuboid(mass: Real, wx: Real, wy: Real, wz: Real) -> Self {
        let c = mass / 12.0;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(
                c * (wy * wy + wz * wz),
                0.0,
                0.0,
                0.0,
                c * (wx * wx + wz * wz),
                0.0,
                0.0,
                0.0,
                c * (wx * wx + wy * wy),
            ),
        }
    }
    /// Mass properties for a uniform-density cylinder (axis along Y).
    pub fn cylinder(mass: Real, radius: Real, height: Real) -> Self {
        let ix = mass * (3.0 * radius * radius + height * height) / 12.0;
        let iy = 0.5 * mass * radius * radius;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ix),
        }
    }
    /// Check if this represents a static (infinite mass) body.
    pub fn is_static(&self) -> bool {
        self.mass <= 0.0
    }
    /// Combine two mass properties (additive).
    pub fn combine(&self, other: &MassProperties) -> MassProperties {
        MassProperties {
            mass: self.mass + other.mass,
            center_of_mass: (self.center_of_mass * self.mass + other.center_of_mass * other.mass)
                / (self.mass + other.mass),
            local_inertia: self.local_inertia + other.local_inertia,
        }
    }
}
impl MassProperties {
    /// Mass properties for a uniform-density cone (apex at top, base at bottom,
    /// axis along +Y). Height = `height`, base radius = `radius`.
    pub fn cone(mass: Real, radius: Real, height: Real) -> Self {
        let ix = mass * (3.0 * radius * radius / 20.0 + height * height * 3.0 / 80.0);
        let iy = 3.0 * mass * radius * radius / 10.0;
        Self {
            mass,
            center_of_mass: Vec3::new(0.0, height * 0.25, 0.0),
            local_inertia: Mat3::new(ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ix),
        }
    }
    /// Mass properties for a hollow spherical shell.
    pub fn spherical_shell(mass: Real, radius: Real) -> Self {
        let i_scalar = (2.0 / 3.0) * mass * radius * radius;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(i_scalar, 0.0, 0.0, 0.0, i_scalar, 0.0, 0.0, 0.0, i_scalar),
        }
    }
    /// Mass properties for a uniform-density capsule (cylinder + 2 hemispheres,
    /// axis along Y).
    pub fn capsule(mass: Real, radius: Real, half_height: Real) -> Self {
        let h = 2.0 * half_height;
        let vol_cyl = std::f64::consts::PI * radius * radius * h;
        let vol_hemi = (2.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        let vol_total = vol_cyl + 2.0 * vol_hemi;
        let density = mass / vol_total;
        let m_cyl = density * vol_cyl;
        let m_hemi = density * vol_hemi;
        let ix_cyl = m_cyl * (3.0 * radius * radius + h * h) / 12.0;
        let iy_cyl = 0.5 * m_cyl * radius * radius;
        let ix_hemi = 83.0 * m_hemi * radius * radius / 320.0;
        let iy_hemi = 2.0 * m_hemi * radius * radius / 5.0;
        let d = half_height + 3.0 * radius / 8.0;
        let ix_hemi_total = ix_hemi + m_hemi * d * d;
        let ix = ix_cyl + 2.0 * ix_hemi_total;
        let iy = iy_cyl + 2.0 * iy_hemi;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ix),
        }
    }
    /// Parallel-axis theorem: shift the inertia tensor by offset `d` from the
    /// center of mass.
    pub fn shifted_inertia(&self, d: Vec3) -> Mat3 {
        let m = self.mass;
        let d2 = d.dot(&d);
        let outer = Mat3::new(
            d.x * d.x,
            d.x * d.y,
            d.x * d.z,
            d.y * d.x,
            d.y * d.y,
            d.y * d.z,
            d.z * d.x,
            d.z * d.y,
            d.z * d.z,
        );
        self.local_inertia + m * (Mat3::identity() * d2 - outer)
    }
}
/// Handle to a rigid body, with generation counter for safe reuse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyHandle {
    /// Index into the body storage.
    pub index: u32,
    /// Generation counter for detecting stale handles.
    pub generation: u32,
}
impl BodyHandle {
    /// Create a new body handle.
    pub fn new(index: u32, generation: u32) -> Self {
        Self { index, generation }
    }
}
impl BodyHandle {
    /// Check if this is a valid (non-null) handle.
    pub fn is_valid(&self) -> bool {
        self.generation > 0 || self.index > 0
    }
    /// Create a null handle.
    pub fn null() -> Self {
        Self {
            index: 0,
            generation: 0,
        }
    }
}
