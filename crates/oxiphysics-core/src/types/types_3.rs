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
    ///
    /// Both parts must be expressed in the same local frame. The result is the
    /// rigid composite of the two parts:
    ///
    /// * the masses add,
    /// * the centre of mass is the mass-weighted mean of the two centres, and
    /// * the inertia tensor about the *combined* centre of mass is the sum of
    ///   each part's own inertia plus its parallel-axis term
    ///   ([`MassProperties::shifted_inertia`]) for the offset of that part's
    ///   centre of mass from the combined one.
    ///
    /// If the total mass is not positive (for example two massless parts) no
    /// mass-weighted centre exists. The result is then a massless property set
    /// located midway between the two centres with a zero inertia tensor, so
    /// the centre of mass is always finite for finite inputs.
    pub fn combine(&self, other: &MassProperties) -> MassProperties {
        let total_mass = self.mass + other.mass;
        if total_mass > 0.0 {
            let center_of_mass =
                (self.center_of_mass * self.mass + other.center_of_mass * other.mass) / total_mass;
            let local_inertia = self.shifted_inertia(self.center_of_mass - center_of_mass)
                + other.shifted_inertia(other.center_of_mass - center_of_mass);
            MassProperties {
                mass: total_mass,
                center_of_mass,
                local_inertia,
            }
        } else {
            MassProperties {
                mass: 0.0,
                center_of_mass: (self.center_of_mass + other.center_of_mass) * 0.5,
                local_inertia: Mat3::zeros(),
            }
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
    /// The null handle: a sentinel that refers to no body.
    ///
    /// Slot arenas hand out generation `0` for the first occupant of a slot
    /// and start at index `0`, so `(0, 0)` is a perfectly ordinary handle (the
    /// very first body of a world). The sentinel therefore uses the largest
    /// index and generation, a pair no arena can realistically issue.
    pub const NULL: Self = Self {
        index: u32::MAX,
        generation: u32::MAX,
    };
    /// Create a new body handle.
    pub fn new(index: u32, generation: u32) -> Self {
        Self { index, generation }
    }
}
impl BodyHandle {
    /// Check if this is a valid (non-null) handle.
    ///
    /// Only [`BodyHandle::NULL`] is invalid. A valid handle may still be stale;
    /// ask the owning arena whether the body is alive.
    pub fn is_valid(&self) -> bool {
        *self != Self::NULL
    }
    /// Create a null handle (see [`BodyHandle::NULL`]).
    pub fn null() -> Self {
        Self::NULL
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ColliderHandle;
    use crate::world::{Body, PhysicsWorld};

    const TOL: Real = 1e-12;

    fn assert_mat_close(actual: &Mat3, expected: &Mat3) {
        for row in 0..3 {
            for col in 0..3 {
                assert!(
                    (actual[(row, col)] - expected[(row, col)]).abs() < TOL,
                    "entry ({row}, {col}): {} vs {}",
                    actual[(row, col)],
                    expected[(row, col)]
                );
            }
        }
    }

    fn assert_vec_close(actual: &Vec3, expected: &Vec3) {
        assert!(
            (actual - expected).norm() < TOL,
            "{actual:?} differs from {expected:?}"
        );
    }

    /// Inertia tensor of point masses about `origin`: `sum m (|r|^2 I - r r^T)`.
    fn point_cloud_inertia(points: &[(Real, Vec3)], origin: &Vec3) -> Mat3 {
        points.iter().fold(Mat3::zeros(), |acc, (mass, position)| {
            let r = position - origin;
            acc + (Mat3::identity() * r.dot(&r) - r * r.transpose()) * *mass
        })
    }

    #[test]
    fn combine_applies_the_parallel_axis_theorem() {
        let left = MassProperties::new(1.0, Vec3::new(-1.0, 0.0, 0.0), Mat3::zeros());
        let right = MassProperties::new(1.0, Vec3::new(1.0, 0.0, 0.0), Mat3::zeros());
        let both = left.combine(&right);
        assert!((both.mass - 2.0).abs() < TOL);
        assert_vec_close(&both.center_of_mass, &Vec3::zeros());
        // Two unit masses on the x axis: nothing about x, 1 + 1 about y and z.
        assert_mat_close(
            &both.local_inertia,
            &Mat3::from_diagonal(&Vec3::new(0.0, 2.0, 2.0)),
        );
    }

    #[test]
    fn combine_matches_the_point_mass_definition_of_inertia() {
        let a = MassProperties::new(2.0, Vec3::new(1.0, 2.0, -0.5), Mat3::zeros());
        let b = MassProperties::new(3.0, Vec3::new(-2.0, 0.5, 1.0), Mat3::zeros());
        let combined = a.combine(&b);
        let expected_com = (a.center_of_mass * 2.0 + b.center_of_mass * 3.0) / 5.0;
        assert_vec_close(&combined.center_of_mass, &expected_com);
        let expected = point_cloud_inertia(
            &[(2.0, a.center_of_mass), (3.0, b.center_of_mass)],
            &expected_com,
        );
        assert_mat_close(&combined.local_inertia, &expected);
    }

    #[test]
    fn combine_keeps_each_part_own_inertia() {
        // Two solid spheres stacked along y: the shift only adds m d^2 about
        // the two axes perpendicular to the stacking direction.
        let sphere = MassProperties::sphere(1.0, 0.5);
        let lower = MassProperties::new(1.0, Vec3::new(0.0, -1.0, 0.0), sphere.local_inertia);
        let upper = MassProperties::new(1.0, Vec3::new(0.0, 1.0, 0.0), sphere.local_inertia);
        let stack = lower.combine(&upper);
        let own = sphere.local_inertia[(0, 0)];
        assert_vec_close(&stack.center_of_mass, &Vec3::zeros());
        assert_mat_close(
            &stack.local_inertia,
            &Mat3::from_diagonal(&Vec3::new(2.0 * own + 2.0, 2.0 * own, 2.0 * own + 2.0)),
        );
    }

    #[test]
    fn combine_is_commutative_and_associative() {
        let a = MassProperties::new(1.5, Vec3::new(0.2, -0.1, 0.4), Mat3::identity() * 0.3);
        let b = MassProperties::new(0.5, Vec3::new(-1.0, 0.7, 0.0), Mat3::identity() * 0.1);
        let c = MassProperties::new(4.0, Vec3::new(0.3, 0.3, -2.0), Mat3::identity() * 0.8);
        let ab = a.combine(&b);
        let ba = b.combine(&a);
        assert!((ab.mass - ba.mass).abs() < TOL);
        assert_vec_close(&ab.center_of_mass, &ba.center_of_mass);
        assert_mat_close(&ab.local_inertia, &ba.local_inertia);
        let left = ab.combine(&c);
        let right = a.combine(&b.combine(&c));
        assert!((left.mass - 6.0).abs() < TOL);
        assert_vec_close(&left.center_of_mass, &right.center_of_mass);
        assert_mat_close(&left.local_inertia, &right.local_inertia);
    }

    #[test]
    fn combining_with_a_massless_part_changes_nothing() {
        let part = MassProperties::new(
            2.0,
            Vec3::new(1.0, -1.0, 0.5),
            Mat3::from_diagonal(&Vec3::new(0.4, 0.5, 0.6)),
        );
        let combined = part.combine(&MassProperties::point_mass(0.0));
        assert!((combined.mass - 2.0).abs() < TOL);
        assert_vec_close(&combined.center_of_mass, &part.center_of_mass);
        assert_mat_close(&combined.local_inertia, &part.local_inertia);
    }

    #[test]
    fn combining_two_massless_parts_has_a_finite_centre() {
        let a = MassProperties::new(0.0, Vec3::new(2.0, 0.0, 0.0), Mat3::zeros());
        let b = MassProperties::new(0.0, Vec3::new(4.0, 2.0, 0.0), Mat3::zeros());
        let combined = a.combine(&b);
        assert_eq!(combined.mass, 0.0);
        assert!(combined.is_static());
        assert!(combined.center_of_mass.iter().all(|c| c.is_finite()));
        assert_vec_close(&combined.center_of_mass, &Vec3::new(3.0, 1.0, 0.0));
        assert_mat_close(&combined.local_inertia, &Mat3::zeros());
    }

    #[test]
    fn null_is_the_only_invalid_body_handle() {
        assert!(!BodyHandle::null().is_valid());
        assert_eq!(BodyHandle::null(), BodyHandle::NULL);
        assert!(!BodyHandle::new(u32::MAX, u32::MAX).is_valid());
        for (index, generation) in [(0, 0), (0, 1), (1, 0), (u32::MAX, 0), (0, u32::MAX)] {
            assert!(
                BodyHandle::new(index, generation).is_valid(),
                "({index}, {generation}) must be a valid handle"
            );
        }
    }

    #[test]
    fn null_is_the_only_invalid_collider_handle() {
        assert!(!ColliderHandle::null().is_valid());
        assert_eq!(ColliderHandle::null(), ColliderHandle::NULL);
        assert!(!ColliderHandle::new(u32::MAX, u32::MAX).is_valid());
        for (index, generation) in [(0, 0), (0, 1), (1, 0), (u32::MAX, 0), (0, u32::MAX)] {
            assert!(
                ColliderHandle::new(index, generation).is_valid(),
                "({index}, {generation}) must be a valid handle"
            );
        }
    }

    #[test]
    fn every_handle_issued_by_the_world_arena_is_valid() {
        let mut world = PhysicsWorld::new();
        let first = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        assert_eq!((first.index, first.generation), (0, 0));
        assert!(first.is_valid(), "the very first body must be addressable");
        assert_ne!(first, BodyHandle::null());
        assert!(world.body(first).is_some());
        assert!(world.body(BodyHandle::null()).is_none());

        let second = world.add_body(Body::new_dynamic([1.0, 0.0, 0.0], 1.0));
        assert!(second.is_valid());
        assert!(world.remove_body(first).is_some());
        let recycled = world.add_body(Body::new_dynamic([2.0, 0.0, 0.0], 1.0));
        assert_eq!(recycled.index, first.index);
        assert_ne!(recycled.generation, first.generation);
        assert!(recycled.is_valid());
        assert!(world.body(first).is_none(), "the old handle is stale");
    }
}
