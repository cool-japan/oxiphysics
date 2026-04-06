//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

/// Collision constraint that prevents a particle from penetrating a plane.
#[derive(Debug, Clone)]
pub struct CollisionConstraint {
    /// Index of the particle.
    pub particle_index: usize,
    /// Plane normal (pointing away from the solid side).
    pub normal: Vec3,
    /// A point on the plane.
    pub point_on_plane: Vec3,
}
impl CollisionConstraint {
    /// Create a new plane-collision constraint.
    pub fn new(particle_index: usize, normal: Vec3, point_on_plane: Vec3) -> Self {
        Self {
            particle_index,
            normal,
            point_on_plane,
        }
    }
}
/// Long Range Attachment constraint (Kim et al. 2012).
///
/// Prevents cloth from stretching beyond a maximum length from a fixed anchor,
/// while allowing free motion within the rest length. Only corrects when the
/// particle is too far from the anchor.
#[derive(Debug, Clone)]
pub struct LraConstraint {
    /// Index of the particle to constrain.
    pub particle: usize,
    /// World-space anchor position.
    pub anchor: Vec3,
    /// Maximum allowed distance from the anchor.
    pub max_distance: Real,
    /// Compliance (XPBD).
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl LraConstraint {
    /// Create a new LRA constraint.
    pub fn new(particle: usize, anchor: Vec3, max_distance: Real, compliance: Real) -> Self {
        Self {
            particle,
            anchor,
            max_distance,
            compliance,
            lambda: 0.0,
        }
    }
    /// Build from current particle position (uses current distance as max).
    pub fn from_particle(
        particle: usize,
        particles: &[SoftParticle],
        anchor: Vec3,
        compliance: Real,
    ) -> Self {
        let dist = (particles[particle].position - anchor).norm();
        Self::new(particle, anchor, dist, compliance)
    }
    /// Reset Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Evaluate constraint violation (positive when too far).
    pub fn constraint_value(&self, particles: &[SoftParticle]) -> Real {
        let dist = (particles[self.particle].position - self.anchor).norm();
        (dist - self.max_distance).max(0.0)
    }
}
/// Surface tension constraint for cloth simulation.
///
/// Minimises local area to simulate surface tension effects. Models cloth
/// wrinkling resistance and thin-film surface energy.
#[derive(Debug, Clone)]
pub struct SurfaceTensionConstraint {
    /// Triangle vertex indices.
    pub indices: [usize; 3],
    /// Rest area of this triangle.
    pub rest_area: Real,
    /// Surface tension coefficient γ (N/m).
    pub gamma: Real,
    /// Compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl SurfaceTensionConstraint {
    /// Create from current particle positions.
    pub fn from_particles(
        indices: [usize; 3],
        particles: &[SoftParticle],
        gamma: Real,
        compliance: Real,
    ) -> Self {
        let p0 = particles[indices[0]].position;
        let p1 = particles[indices[1]].position;
        let p2 = particles[indices[2]].position;
        let rest_area = 0.5 * (p1 - p0).cross(&(p2 - p0)).norm();
        Self {
            indices,
            rest_area,
            gamma,
            compliance,
            lambda: 0.0,
        }
    }
    /// Reset Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Current triangle area.
    pub fn current_area(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> Real {
        0.5 * (p1 - p0).cross(&(p2 - p0)).norm()
    }
    /// Surface energy: E = γ * (A - A0).
    pub fn surface_energy(&self, particles: &[SoftParticle]) -> Real {
        let [i0, i1, i2] = self.indices;
        let area = Self::current_area(
            &particles[i0].position,
            &particles[i1].position,
            &particles[i2].position,
        );
        self.gamma * (area - self.rest_area)
    }
}
/// Shape matching constraint (Müller et al. 2005).
///
/// Preserves the shape of a group of particles by computing the optimal
/// rotation from the rest shape to the current configuration.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ShapeMatchingConstraint {
    /// Indices of the particles in this shape.
    pub indices: Vec<usize>,
    /// Rest positions relative to center of mass.
    pub rest_relative: Vec<Vec3>,
    /// Rest center of mass.
    pub rest_com: Vec3,
    /// Stiffness (0-1, where 1 is rigid).
    pub stiffness: Real,
}
impl ShapeMatchingConstraint {
    /// Create a shape matching constraint from the current configuration.
    #[allow(dead_code)]
    pub fn from_particles(
        indices: Vec<usize>,
        particles: &[SoftParticle],
        stiffness: Real,
    ) -> Self {
        let n = indices.len();
        if n == 0 {
            return Self {
                indices,
                rest_relative: Vec::new(),
                rest_com: Vec3::zeros(),
                stiffness,
            };
        }
        let total_inv_mass: Real = indices
            .iter()
            .map(|&i| {
                if particles[i].inverse_mass > 0.0 {
                    1.0 / particles[i].inverse_mass
                } else {
                    1e6
                }
            })
            .sum();
        let mut com = Vec3::zeros();
        for &i in &indices {
            let mass = if particles[i].inverse_mass > 0.0 {
                1.0 / particles[i].inverse_mass
            } else {
                1e6
            };
            com += particles[i].position * mass;
        }
        com /= total_inv_mass;
        let rest_relative: Vec<Vec3> = indices
            .iter()
            .map(|&i| particles[i].position - com)
            .collect();
        Self {
            indices,
            rest_relative,
            rest_com: com,
            stiffness,
        }
    }
    /// Compute current center of mass.
    pub(super) fn current_com(&self, particles: &[SoftParticle]) -> Vec3 {
        let n = self.indices.len();
        if n == 0 {
            return Vec3::zeros();
        }
        let mut total_mass = 0.0_f64;
        let mut com = Vec3::zeros();
        for &i in &self.indices {
            let mass = if particles[i].inverse_mass > 0.0 {
                1.0 / particles[i].inverse_mass
            } else {
                1e6
            };
            com += particles[i].position * mass;
            total_mass += mass;
        }
        if total_mass > 1e-14 {
            com / total_mass
        } else {
            Vec3::zeros()
        }
    }
    /// Compute the optimal rotation using the polar decomposition
    /// via cross-covariance matrix.
    ///
    /// Returns the rotation matrix as a 3x3 array (row-major).
    pub(super) fn compute_rotation(&self, particles: &[SoftParticle]) -> [Real; 9] {
        let com = self.current_com(particles);
        let mut a = [0.0_f64; 9];
        for (k, &idx) in self.indices.iter().enumerate() {
            let mass = if particles[idx].inverse_mass > 0.0 {
                1.0 / particles[idx].inverse_mass
            } else {
                1e6
            };
            let p = particles[idx].position - com;
            let q = &self.rest_relative[k];
            a[0] += mass * p.x * q.x;
            a[1] += mass * p.x * q.y;
            a[2] += mass * p.x * q.z;
            a[3] += mass * p.y * q.x;
            a[4] += mass * p.y * q.y;
            a[5] += mass * p.y * q.z;
            a[6] += mass * p.z * q.x;
            a[7] += mass * p.z * q.y;
            a[8] += mass * p.z * q.z;
        }
        let col0 = Vec3::new(a[0], a[3], a[6]);
        let col1_raw = Vec3::new(a[1], a[4], a[7]);
        let col2_raw = Vec3::new(a[2], a[5], a[8]);
        let len0 = col0.norm();
        if len0 < 1e-14 {
            return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        }
        let e0 = col0 / len0;
        let proj1 = col1_raw - e0 * e0.dot(&col1_raw);
        let len1 = proj1.norm();
        if len1 < 1e-14 {
            return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        }
        let e1 = proj1 / len1;
        let e2_raw = col2_raw - e0 * e0.dot(&col2_raw) - e1 * e1.dot(&col2_raw);
        let len2 = e2_raw.norm();
        if len2 < 1e-14 {
            return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        }
        let e2 = e2_raw / len2;
        let det_sign = e0.dot(&e1.cross(&e2));
        let e2 = if det_sign < 0.0 { -e2 } else { e2 };
        [e0.x, e1.x, e2.x, e0.y, e1.y, e2.y, e0.z, e1.z, e2.z]
    }
}
/// Strain limiting constraint for a triangle element.
///
/// Limits the maximum stretch and compression of a triangle to prevent
/// excessive deformation. Based on principal strain analysis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct StrainLimitingConstraint {
    /// Indices of the three vertices.
    pub indices: [usize; 3],
    /// Maximum allowed stretch ratio (e.g. 1.1 for 10% stretch).
    pub max_stretch: Real,
    /// Minimum allowed compression ratio (e.g. 0.9 for 10% compression).
    pub min_compression: Real,
    /// Inverse rest transform (2x2 stored as \[a, b, c, d\]).
    pub(super) inv_rest: [Real; 4],
    /// Compliance (XPBD).
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl StrainLimitingConstraint {
    /// Create a strain limiting constraint from the current configuration.
    #[allow(dead_code)]
    pub fn from_particles(
        indices: [usize; 3],
        particles: &[SoftParticle],
        max_stretch: Real,
        min_compression: Real,
        compliance: Real,
    ) -> Self {
        let p0 = particles[indices[0]].position;
        let p1 = particles[indices[1]].position;
        let p2 = particles[indices[2]].position;
        let inv_rest = Self::compute_inv_rest(&p0, &p1, &p2);
        Self {
            indices,
            max_stretch,
            min_compression,
            inv_rest,
            compliance,
            lambda: 0.0,
        }
    }
    /// Reset the Lagrange multiplier.
    #[allow(dead_code)]
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Compute the inverse of the 2D rest matrix.
    fn compute_inv_rest(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> [Real; 4] {
        let d1 = p1 - p0;
        let d2 = p2 - p0;
        let e1 = d1;
        let e1_len = e1.norm();
        if e1_len < 1e-14 {
            return [1.0, 0.0, 0.0, 1.0];
        }
        let u = e1 / e1_len;
        let normal = d1.cross(&d2);
        let n_len = normal.norm();
        if n_len < 1e-14 {
            return [1.0, 0.0, 0.0, 1.0];
        }
        let v = normal.cross(&d1).normalize();
        let a = d1.dot(&u);
        let b = d1.dot(&v);
        let c = d2.dot(&u);
        let d = d2.dot(&v);
        let det = a * d - b * c;
        if det.abs() < 1e-14 {
            return [1.0, 0.0, 0.0, 1.0];
        }
        [d / det, -c / det, -b / det, a / det]
    }
}
/// Neo-Hookean strain energy constraint for a tetrahedron (FEM-based XPBD).
///
/// Uses the linearised Neo-Hookean model:
///   Ψ = (mu/2)(I1 - 3) + (lambda/2)(J-1)^2
/// where I1 = tr(F^T F), J = det(F), and F is the deformation gradient.
///
/// Indices: \[i0, i1, i2, i3\] — tet vertices.
#[derive(Debug, Clone)]
pub struct NeoHookeanConstraint {
    /// Indices of the four tet vertices.
    pub indices: [usize; 4],
    /// Shear modulus (mu).
    pub mu: Real,
    /// First Lamé parameter (lambda).
    pub lambda_lame: Real,
    /// Compliance for volume part.
    pub compliance: Real,
    /// Inverse rest-shape matrix Dm^-1 (stored as 9 entries, column-major 3x3).
    pub(super) dm_inv: [Real; 9],
    /// Rest volume of tetrahedron.
    pub rest_volume: Real,
    /// Accumulated Lagrange multiplier (deviatoric).
    pub(super) lambda_dev: Real,
    /// Accumulated Lagrange multiplier (hydrostatic).
    pub(super) lambda_vol: Real,
}
#[allow(clippy::too_many_arguments)]
impl NeoHookeanConstraint {
    /// Build constraint from current particle positions.
    pub fn from_particles(
        indices: [usize; 4],
        particles: &[SoftParticle],
        mu: Real,
        lambda_lame: Real,
        compliance: Real,
    ) -> Self {
        let p0 = particles[indices[0]].position;
        let p1 = particles[indices[1]].position;
        let p2 = particles[indices[2]].position;
        let p3 = particles[indices[3]].position;
        let dm = Self::compute_shape_matrix(&p0, &p1, &p2, &p3);
        let dm_inv = mat3_inverse(dm).unwrap_or([0.0; 9]);
        let rest_volume = VolumeConstraint::compute_tet_volume(&p0, &p1, &p2, &p3).abs();
        Self {
            indices,
            mu,
            lambda_lame,
            compliance,
            dm_inv,
            rest_volume,
            lambda_dev: 0.0,
            lambda_vol: 0.0,
        }
    }
    /// Reset Lagrange multipliers.
    pub fn reset_lambda(&mut self) {
        self.lambda_dev = 0.0;
        self.lambda_vol = 0.0;
    }
    /// Compute deformation gradient F = Ds * Dm_inv.
    pub fn deformation_gradient(&self, particles: &[SoftParticle]) -> [Real; 9] {
        let [i0, i1, i2, i3] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let p3 = particles[i3].position;
        let ds = Self::compute_shape_matrix(&p0, &p1, &p2, &p3);
        mat3_mul(ds, self.dm_inv)
    }
    /// Compute I1 = tr(F^T F).
    pub fn compute_i1(f: [Real; 9]) -> Real {
        f.iter().map(|&x| x * x).sum()
    }
    /// Compute J = det(F).
    pub fn compute_j(f: [Real; 9]) -> Real {
        mat3_det(f)
    }
    fn compute_shape_matrix(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> [Real; 9] {
        let c0 = *p1 - *p0;
        let c1 = *p2 - *p0;
        let c2 = *p3 - *p0;
        [c0.x, c0.y, c0.z, c1.x, c1.y, c1.z, c2.x, c2.y, c2.z]
    }
}
/// Area conservation constraint for a single triangle.
///
/// Preserves the area of a triangle formed by three particles.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AreaConstraint {
    /// Indices of the three vertices.
    pub indices: [usize; 3],
    /// Rest area.
    pub rest_area: Real,
    /// Compliance (XPBD).
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl AreaConstraint {
    /// Create an area constraint with an explicit rest area.
    #[allow(dead_code)]
    pub fn new(indices: [usize; 3], rest_area: Real, compliance: Real) -> Self {
        Self {
            indices,
            rest_area,
            compliance,
            lambda: 0.0,
        }
    }
    /// Build an area constraint from the current configuration.
    #[allow(dead_code)]
    pub fn from_particles(
        indices: [usize; 3],
        particles: &[SoftParticle],
        compliance: Real,
    ) -> Self {
        let area = Self::compute_triangle_area(
            &particles[indices[0]].position,
            &particles[indices[1]].position,
            &particles[indices[2]].position,
        );
        Self::new(indices, area, compliance)
    }
    /// Reset the Lagrange multiplier.
    #[allow(dead_code)]
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Compute the area of a triangle.
    pub(super) fn compute_triangle_area(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> Real {
        let e1 = p1 - p0;
        let e2 = p2 - p0;
        0.5 * e1.cross(&e2).norm()
    }
}
/// Binding between a soft-body particle and a rigid-body anchor point.
///
/// Pulls the particle towards a world-space point that may be updated each
/// frame (e.g. driven by a rigid body transform).
#[derive(Debug, Clone)]
pub struct RigidBindingConstraint {
    /// Particle index.
    pub particle: usize,
    /// World-space target position (updated externally).
    pub target: Vec3,
    /// Compliance (0 = hard attachment).
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl RigidBindingConstraint {
    /// Create a new rigid binding constraint.
    pub fn new(particle: usize, target: Vec3, compliance: Real) -> Self {
        Self {
            particle,
            target,
            compliance,
            lambda: 0.0,
        }
    }
    /// Reset Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Update the target position (call before each solver iteration).
    pub fn set_target(&mut self, target: Vec3) {
        self.target = target;
    }
    /// Residual (distance to target).
    pub fn residual(&self, particles: &[SoftParticle]) -> Real {
        (particles[self.particle].position - self.target).norm()
    }
}
/// Motor constraint for soft actuators.
///
/// Drives a pair of particles towards a target length (like a muscle fibre).
/// Uses XPBD with a target rest length that can be set dynamically.
#[derive(Debug, Clone)]
pub struct MotorConstraint {
    /// First particle index.
    pub i: usize,
    /// Second particle index.
    pub j: usize,
    /// Current target length (driven by actuation signal).
    pub target_length: Real,
    /// Natural rest length (fully relaxed).
    pub rest_length: Real,
    /// Maximum contraction ratio (e.g. 0.6 means 60% of rest).
    pub min_ratio: Real,
    /// Compliance.
    pub compliance: Real,
    /// Current actuation in \[0, 1\]: 0 = relaxed, 1 = fully contracted.
    pub actuation: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl MotorConstraint {
    /// Create a new motor constraint.
    pub fn new(i: usize, j: usize, rest_length: Real, min_ratio: Real, compliance: Real) -> Self {
        Self {
            i,
            j,
            target_length: rest_length,
            rest_length,
            min_ratio,
            compliance,
            actuation: 0.0,
            lambda: 0.0,
        }
    }
    /// Build from current particle distance.
    pub fn from_particles(
        i: usize,
        j: usize,
        particles: &[SoftParticle],
        min_ratio: Real,
        compliance: Real,
    ) -> Self {
        let rest = (particles[i].position - particles[j].position).norm();
        Self::new(i, j, rest, min_ratio, compliance)
    }
    /// Reset Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Set the actuation signal \[0, 1\] and update target length.
    pub fn set_actuation(&mut self, actuation: Real) {
        self.actuation = actuation.clamp(0.0, 1.0);
        let min_len = self.rest_length * self.min_ratio;
        self.target_length = self.rest_length - (self.rest_length - min_len) * self.actuation;
    }
    /// Current contraction ratio.
    pub fn contraction_ratio(&self) -> Real {
        if self.rest_length < 1e-12 {
            return 1.0;
        }
        self.target_length / self.rest_length
    }
}
/// Volume-preservation constraint for a tetrahedral mesh.
///
/// Operates over all tetrahedra whose indices are stored in the constraint.
#[derive(Debug, Clone)]
pub struct VolumeConstraint {
    /// Indices of the four vertices forming the tetrahedron.
    pub indices: [usize; 4],
    /// Rest volume (positive).
    pub rest_volume: Real,
    /// Compliance (XPBD).
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl VolumeConstraint {
    /// Create a volume constraint with an explicit rest volume.
    pub fn new(indices: [usize; 4], rest_volume: Real, compliance: Real) -> Self {
        Self {
            indices,
            rest_volume,
            compliance,
            lambda: 0.0,
        }
    }
    /// Build a volume constraint whose rest volume matches the current
    /// configuration.
    pub fn from_particles(
        indices: [usize; 4],
        particles: &[SoftParticle],
        compliance: Real,
    ) -> Self {
        let vol = Self::compute_tet_volume(
            &particles[indices[0]].position,
            &particles[indices[1]].position,
            &particles[indices[2]].position,
            &particles[indices[3]].position,
        );
        Self::new(indices, vol, compliance)
    }
    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Signed volume of a tetrahedron.
    pub fn compute_tet_volume(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> Real {
        let a = p1 - p0;
        let b = p2 - p0;
        let c = p3 - p0;
        a.dot(&b.cross(&c)) / 6.0
    }
}
/// XPBD distance constraint between two particles.
///
/// Maintains a rest length with optional compliance (softness).
#[derive(Debug, Clone)]
pub struct DistanceConstraint {
    /// Index of the first particle.
    pub i: usize,
    /// Index of the second particle.
    pub j: usize,
    /// Rest length.
    pub rest_length: Real,
    /// Compliance (inverse stiffness). 0 = perfectly rigid.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier (reset each sub-step sequence).
    pub(super) lambda: Real,
}
impl DistanceConstraint {
    /// Create a new distance constraint.
    pub fn new(i: usize, j: usize, rest_length: Real, compliance: Real) -> Self {
        Self {
            i,
            j,
            rest_length,
            compliance,
            lambda: 0.0,
        }
    }
    /// Build a constraint whose rest length is the current distance between the
    /// two particles.
    pub fn from_particles(
        i: usize,
        j: usize,
        particles: &[SoftParticle],
        compliance: Real,
    ) -> Self {
        let rest = (particles[i].position - particles[j].position).norm();
        Self::new(i, j, rest, compliance)
    }
    /// Reset the Lagrange multiplier (call before a new solver iteration
    /// sequence).
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
}
/// Isometric bending constraint (Bergou et al. 2006).
///
/// Uses a quadratic energy based on the Laplacian of the surface positions,
/// which is more robust than dihedral-angle bending for small angles.
///
/// For two triangles sharing edge (p0,p1) with opposite vertices p2,p3:
/// `C = (p0*Q[0] + p1*Q[1] + p2*Q[2] + p3*Q[3])`
/// where Q are cotangent-weighted coefficients.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct IsometricBendingConstraint {
    /// Indices of the four particles: shared edge (i0, i1), wing vertices (i2, i3).
    pub indices: [usize; 4],
    /// Stiffness coefficient for bending.
    pub stiffness: Real,
    /// Compliance (XPBD).
    pub compliance: Real,
    /// Rest state matrix Q (4x4 symmetric, stored as flat array).
    pub(super) q_matrix: [Real; 16],
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl IsometricBendingConstraint {
    /// Create an isometric bending constraint from the current configuration.
    #[allow(dead_code)]
    pub fn from_particles(
        indices: [usize; 4],
        particles: &[SoftParticle],
        stiffness: Real,
        compliance: Real,
    ) -> Self {
        let p0 = particles[indices[0]].position;
        let p1 = particles[indices[1]].position;
        let p2 = particles[indices[2]].position;
        let p3 = particles[indices[3]].position;
        let q_matrix = Self::compute_q_matrix(&p0, &p1, &p2, &p3);
        Self {
            indices,
            stiffness,
            compliance,
            q_matrix,
            lambda: 0.0,
        }
    }
    /// Reset the Lagrange multiplier.
    #[allow(dead_code)]
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Compute the Q matrix from cotangent weights.
    ///
    /// Uses the Bergou et al. formulation where the constraint is based on
    /// the discrete Laplacian. The key idea is that the constraint vector
    /// K = \[c03+c06, c02+c05, -c01, -c04\] satisfies K.p = 0 at rest
    /// (when the surface is flat), and Q = K^T K / (4A) encodes the energy.
    fn compute_q_matrix(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> [Real; 16] {
        let e = p1 - p0;
        let e0_2 = p2 - p0;
        let e1_2 = p2 - p1;
        let e0_3 = p3 - p0;
        let e1_3 = p3 - p1;
        let cot = |a: &Vec3, b: &Vec3| -> Real {
            let cross_norm = a.cross(b).norm();
            if cross_norm < 1e-14 {
                return 0.0;
            }
            a.dot(b) / cross_norm
        };
        let c01 = cot(&e0_2, &e1_2);
        let c02 = cot(&e, &e0_2);
        let c03 = cot(&(-e), &e1_2);
        let c04 = cot(&e0_3, &e1_3);
        let c05 = cot(&e, &e0_3);
        let c06 = cot(&(-e), &e1_3);
        let a1 = 0.5 * e0_2.cross(&e).norm();
        let a2 = 0.5 * e0_3.cross(&e).norm();
        let a_total = a1 + a2;
        if a_total < 1e-14 {
            return [0.0; 16];
        }
        let kk = [c03 + c06, c02 + c05, -c01, -c04];
        let mut q = [0.0_f64; 16];
        for row in 0..4 {
            for col in 0..4 {
                q[row * 4 + col] = kk[row] * kk[col] / (4.0 * a_total);
            }
        }
        q
    }
}
/// Improved isometric bending constraint that applies gradients to all four
/// vertices (shared edge + two wing vertices), not just the wings.
///
/// Uses the quadratic form from Bergou et al. 2006 for better convergence.
#[derive(Debug, Clone)]
pub struct IsometricBendingConstraintV2 {
    /// \[i0, i1, i2, i3\]: shared edge (i0,i1), wing vertices (i2,i3).
    pub indices: [usize; 4],
    /// Bending stiffness.
    pub stiffness: Real,
    /// Compliance.
    pub compliance: Real,
    /// 4×4 Q matrix (flat, row-major).
    pub(super) q_matrix: [Real; 16],
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl IsometricBendingConstraintV2 {
    /// Build from current positions.
    pub fn from_particles(
        indices: [usize; 4],
        particles: &[SoftParticle],
        stiffness: Real,
        compliance: Real,
    ) -> Self {
        let p = [
            particles[indices[0]].position,
            particles[indices[1]].position,
            particles[indices[2]].position,
            particles[indices[3]].position,
        ];
        let q_matrix = Self::build_q(&p[0], &p[1], &p[2], &p[3]);
        Self {
            indices,
            stiffness,
            compliance,
            q_matrix,
            lambda: 0.0,
        }
    }
    /// Reset Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Quadratic bending energy: E = (stiffness/2) * sum_i p_i^T Q p_i.
    pub fn bending_energy(&self, particles: &[SoftParticle]) -> Real {
        let ps: Vec<Vec3> = self
            .indices
            .iter()
            .map(|&i| particles[i].position)
            .collect();
        let mut energy = 0.0;
        for axis in 0..3 {
            let coords: [Real; 4] = std::array::from_fn(|k| match axis {
                0 => ps[k].x,
                1 => ps[k].y,
                _ => ps[k].z,
            });
            for i in 0..4 {
                for j in 0..4 {
                    energy += self.q_matrix[i * 4 + j] * coords[i] * coords[j];
                }
            }
        }
        0.5 * self.stiffness * energy
    }
    fn cot_angle(a: &Vec3, b: &Vec3) -> Real {
        let cross = a.cross(b).norm();
        if cross < 1e-14 {
            return 0.0;
        }
        a.dot(b) / cross
    }
    fn build_q(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> [Real; 16] {
        let e0 = *p2 - *p0;
        let e1 = *p2 - *p1;
        let e2 = *p3 - *p0;
        let e3 = *p3 - *p1;
        let c01 = Self::cot_angle(&e0, &e1);
        let c02 = Self::cot_angle(&e2, &e3);
        let k = [c01 + c02, -(c01 + c02), c01, -c01 + c02 - c02, c02, -c02];
        let coeff = [k[0], k[1], k[2], k[3]];
        let mut q = [0.0f64; 16];
        for i in 0..4 {
            for j in 0..4 {
                q[i * 4 + j] = coeff[i] * coeff[j];
            }
        }
        q
    }
}
/// Balloon pressure constraint for inflatable objects.
///
/// Enforces a target enclosed volume by applying outward pressure on the
/// surface triangles. Indices specify a surface triangle (i0, i1, i2).
#[derive(Debug, Clone)]
pub struct BalloonPressureConstraint {
    /// Surface triangle indices.
    pub indices: [usize; 3],
    /// Target pressure (multiplied by current volume deviation).
    pub pressure: Real,
    /// Rest volume reference (of the entire balloon).
    pub rest_volume: Real,
    /// Compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl BalloonPressureConstraint {
    /// Create a new balloon pressure constraint.
    pub fn new(indices: [usize; 3], pressure: Real, rest_volume: Real, compliance: Real) -> Self {
        Self {
            indices,
            pressure,
            rest_volume,
            compliance,
            lambda: 0.0,
        }
    }
    /// Reset Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Compute signed triangle area contribution to enclosed volume.
    /// Divergence theorem: V = (1/6) * sum over faces of (p0 · (p1 × p2)).
    pub fn triangle_volume_contribution(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> Real {
        p0.dot(&p1.cross(p2)) / 6.0
    }
    /// Outward normal of the triangle (unnormalized, magnitude = area).
    pub(super) fn triangle_normal_area(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> Vec3 {
        (p1 - p0).cross(&(p2 - p0))
    }
}
/// Dihedral-angle bending constraint for two triangles sharing an edge.
///
/// Vertices (p0, p1) are the shared edge; p2 and p3 are the opposing vertices.
#[derive(Debug, Clone)]
pub struct BendingConstraint {
    /// Indices of the four particles: shared edge (i0, i1), wing vertices (i2, i3).
    pub indices: [usize; 4],
    /// Rest dihedral angle.
    pub rest_angle: Real,
    /// Compliance (XPBD).
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    pub(super) lambda: Real,
}
impl BendingConstraint {
    /// Create a bending constraint with the given rest angle.
    pub fn new(indices: [usize; 4], rest_angle: Real, compliance: Real) -> Self {
        Self {
            indices,
            rest_angle,
            compliance,
            lambda: 0.0,
        }
    }
    /// Build a bending constraint whose rest angle matches the current
    /// configuration.
    pub fn from_particles(
        indices: [usize; 4],
        particles: &[SoftParticle],
        compliance: Real,
    ) -> Self {
        let angle = Self::compute_dihedral(
            &particles[indices[0]].position,
            &particles[indices[1]].position,
            &particles[indices[2]].position,
            &particles[indices[3]].position,
        );
        Self::new(indices, angle, compliance)
    }
    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Compute the dihedral angle between two triangles sharing edge (p0,p1).
    pub(super) fn compute_dihedral(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> Real {
        let e = p1 - p0;
        let n1 = (p2 - p0).cross(&e);
        let n2 = (p3 - p0).cross(&e);
        let n1_len = n1.norm();
        let n2_len = n2.norm();
        if n1_len < 1e-12 || n2_len < 1e-12 {
            return 0.0;
        }
        let n1 = n1 / n1_len;
        let n2 = n2 / n2_len;
        let cos_a = n1.dot(&n2).clamp(-1.0, 1.0);
        cos_a.acos()
    }
}
