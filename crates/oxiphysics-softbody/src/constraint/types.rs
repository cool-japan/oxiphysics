//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::particle::SoftParticle;
use oxiphysics_core::math::{Mat3, Real, Vec3};

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
/// Stand-in mass of a pinned particle (inverse mass zero) in the weighted
/// averages of [`ShapeMatchingConstraint`]: large, so that a pinned particle
/// dominates the centre of mass and the orientation of its group, yet finite.
const PINNED_PARTICLE_MASS: Real = 1.0e6;

/// Shape matching constraint (Müller et al. 2005).
///
/// Preserves the shape of a group of particles by computing the optimal
/// rotation from the rest shape to the current configuration.
///
/// With the offsets `p_i = x_i - c` of the particles from their current centre
/// of mass `c` and the rest offsets `q_i`, the optimal rotation is the rotation
/// factor `R` of the polar decomposition `A = R S` of the mass-weighted
/// cross-covariance matrix `A = sum_i m_i p_i q_i^T`.  It is computed with a
/// singular value decomposition rather than an orthogonalisation of the columns
/// of `A`, so a resting group is not rotated (`A` is then the symmetric rest
/// covariance matrix, whose rotation factor is the identity) and groups whose
/// particles lie in a plane, such as a cloth patch (`A` has rank 2), keep a well
/// defined orientation.  The rotation of a collinear group (rank 1) is only
/// determined up to a twist about its line, which moves no goal position.
/// `project` moves every unpinned particle by `stiffness` times the vector to its
/// goal position `R q_i + c`; without pinned particles the mass-weighted
/// centroid is unchanged.
#[derive(Debug, Clone)]
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
        let total_mass: Real = indices
            .iter()
            .map(|&i| Self::particle_mass(&particles[i]))
            .sum();
        let mut com = Vec3::zeros();
        for &i in &indices {
            com += particles[i].position * Self::particle_mass(&particles[i]);
        }
        com /= total_mass;
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
    /// Mass of a particle as used by the weighted averages of the constraint:
    /// `1 / inverse_mass`, or [`PINNED_PARTICLE_MASS`] for a pinned particle.
    fn particle_mass(particle: &SoftParticle) -> Real {
        if particle.inverse_mass > 0.0 {
            1.0 / particle.inverse_mass
        } else {
            PINNED_PARTICLE_MASS
        }
    }
    /// Compute current center of mass.
    pub(super) fn current_com(&self, particles: &[SoftParticle]) -> Vec3 {
        let mut total_mass = 0.0_f64;
        let mut com = Vec3::zeros();
        for &i in &self.indices {
            let mass = Self::particle_mass(&particles[i]);
            com += particles[i].position * mass;
            total_mass += mass;
        }
        if total_mass > 0.0 {
            com / total_mass
        } else {
            Vec3::zeros()
        }
    }
    /// Rotation `R` of the polar decomposition `A = R S` of the mass-weighted
    /// cross-covariance matrix `A = sum_i m_i (x_i - com) q_i^T`, where `com` is
    /// the current centre of mass (see [`current_com`](Self::current_com)) and
    /// `q_i` the rest offsets.
    ///
    /// `R` is the proper rotation that maps the rest offsets onto the current
    /// ones as closely as a rotation can: `x_i - com ~ R q_i`.  It is the
    /// identity for a group at rest and when the group carries no orientation (a
    /// single particle, coincident particles, non-finite positions); for a
    /// collinear group it is fixed only up to a twist about the line.
    pub(super) fn compute_rotation(&self, particles: &[SoftParticle], com: &Vec3) -> Mat3 {
        let mut covariance = Mat3::zeros();
        for (&idx, q) in self.indices.iter().zip(&self.rest_relative) {
            let particle = &particles[idx];
            let p = particle.position - com;
            covariance += (p * q.transpose()) * Self::particle_mass(particle);
        }
        polar_rotation(&covariance)
    }
}
/// Rotation factor `R` of the polar decomposition `a = R S` (`S` symmetric
/// positive semi-definite): the proper rotation (`R^T R = I`, `det R = +1`) that
/// is closest to `a` in the Frobenius norm.
///
/// With the singular value decomposition `a = U diag(s) V^T` the rotation is
/// `U V^T`.  If that product has determinant `-1` (the matrix mirrors the
/// shape), the closest proper rotation is obtained by negating the column of `U`
/// that belongs to the smallest singular value.  No inverse of `a` is needed, so
/// the factor is also defined for a matrix of rank 2, which is what a planar
/// group has (all `p_i` and `q_i` lie in a plane): the smallest singular value is
/// then zero and the sign choice yields the unique proper rotation that agrees
/// with `a` on the plane, including its orientation when the plane is turned
/// over.
///
/// `a` is normalised by its largest entry first, so the result does not depend
/// on the scale of the masses and extents.  The identity is returned for the
/// zero matrix and for a matrix with a non-finite entry: neither carries an
/// orientation to extract.
pub(super) fn polar_rotation(a: &Mat3) -> Mat3 {
    if a.iter().any(|entry| !entry.is_finite()) {
        return Mat3::identity();
    }
    let scale = a
        .iter()
        .fold(0.0_f64, |largest, entry| largest.max(entry.abs()));
    if scale <= 0.0 {
        return Mat3::identity();
    }
    let svd = (a / scale).svd(true, true);
    let (Some(mut u), Some(v_t)) = (svd.u, svd.v_t) else {
        return Mat3::identity();
    };
    let rotation = u * v_t;
    if rotation.determinant() >= 0.0 {
        return rotation;
    }
    let smallest = svd
        .singular_values
        .iter()
        .enumerate()
        .min_by(|(_, first), (_, second)| first.total_cmp(second))
        .map_or(2, |(index, _)| index);
    u.column_mut(smallest).neg_mut();
    u * v_t
}
/// Strain limiting constraint for a triangle element.
///
/// Limits the principal stretches of a triangle to `[min_compression,
/// max_stretch]` (relative to its rest shape) to prevent excessive deformation.
///
/// In its own plane the triangle has the deformation gradient `F = Ds Dm^-1` with
/// the singular value decomposition `F = U diag(s1, s2) V^T`.  `project` clamps
/// `s1` and `s2` into the limits and moves the vertices to the edge vectors of
/// the clamped gradient `F* = U diag(clamp(s1), clamp(s2)) V^T`, spreading the
/// motion by the inverse masses so that the mass-weighted centroid stays fixed
/// (pinned vertices do not move).  A stretched triangle shrinks only as far as the
/// limit requires, an over-compressed one expands, and a triangle inside its
/// limits is left alone.
///
/// The projection is stateless: there is no Lagrange multiplier to accumulate.  A
/// positive `compliance` softens it in the XPBD manner (see
/// [`SoftConstraint::project`](super::functions::SoftConstraint::project)).
#[derive(Debug, Clone)]
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
}
impl StrainLimitingConstraint {
    /// Create a strain limiting constraint from the current configuration.
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
        }
    }
    /// Reset the accumulated constraint state before a new sub-step.
    ///
    /// The strain limit is a stateless projection, so there is nothing to reset;
    /// the method exists so that solvers can reset every constraint uniformly.
    pub fn reset_lambda(&mut self) {}
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
/// Stable Neo-Hookean strain-energy constraint for a tetrahedron (XPBD).
///
/// Follows the constraint formulation of Macklin and Müller, "A Constraint-based
/// Formulation of Stable Neo-Hookean Materials" (2021).  With the deformation
/// gradient `F = Ds Dm^-1`, `I1 = tr(F^T F)`, `J = det(F)` and the rest volume
/// `V0`, the element energy `V0 [ (mu/2)(I1 - 3) + (K/2)(J - gamma)^2 ]` is split
/// into two XPBD constraints that are projected together:
///
/// * deviatoric: `C_D = sqrt(I1)` with compliance `1 / (mu V0)`,
/// * hydrostatic: `C_H = J - gamma` with `gamma = 1 + mu / K` and compliance
///   `1 / (K V0) + compliance`.
///
/// `K = lambda + mu` is the hydrostatic modulus
/// ([`hydrostatic_modulus`](Self::hydrostatic_modulus)).  Compared with using
/// `lambda` directly, this makes the first Piola stress `P = mu F + K (J - gamma)
/// cof(F)` equal to Hooke's law with the Lamé parameters `(mu, lambda)` for
/// small strains (`P_xx = (lambda + 2 mu) eps`, `P_yy = lambda eps` under
/// uniaxial strain), keeps the rest state stable for every Poisson ratio below
/// 1/2, and stays finite for `lambda = 0`.  With `lambda` itself the effective
/// bulk modulus would be `lambda - mu/3` and the rest state would be unstable
/// for Poisson ratios below 1/8.
///
/// `C_D` equals `sqrt(3)` at rest, so the two constraint forces only cancel
/// when both are solved as one coupled 2 x 2 system; solving them one after the
/// other makes stiff elements drift away from, and eventually collapse towards,
/// zero volume.  The coupled projection leaves the rest state exactly
/// unchanged, conserves the mass-weighted centroid, and pushes inverted
/// elements back to positive volume.
///
/// The XPBD update is a Gauss-Newton step of the element's implicit-Euler
/// problem.  Both constraints are strongly nonlinear in the vertex positions, so
/// when `mu V0 dt_sub^2` is large compared with the element mass (a stiff
/// element, a long sub-step) the full step overshoots and can collapse or invert
/// the element.  `project` therefore halves the step until it sufficiently
/// decreases the element-local incremental potential (Armijo backtracking).  For
/// the usual stiffness and sub-step ranges the full XPBD step passes that test
/// and nothing changes; for stiff elements it makes repeated projections (solver
/// iterations, or `reset_lambda` followed by `project`) relax the element
/// monotonically.  One iteration per sub-step is still only an approximation of
/// the implicit step, so stiff materials call for sub-steps or several
/// iterations.
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
    /// Extra compliance of the hydrostatic (volume) constraint in `1/J`, added in
    /// series to its physical compliance `1 / (K V0)` with `K = lambda + mu`;
    /// `0` gives the pure Neo-Hookean response.
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
        self.deformation_gradient_of(&self.indices.map(|index| particles[index].position))
    }
    /// Deformation gradient F = Ds * Dm_inv of the element whose four corners
    /// are at `positions`, given in the order of `indices`.
    pub(super) fn deformation_gradient_of(&self, positions: &[Vec3; 4]) -> [Real; 9] {
        let ds =
            Self::compute_shape_matrix(&positions[0], &positions[1], &positions[2], &positions[3]);
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
    pub fn new(indices: [usize; 3], rest_area: Real, compliance: Real) -> Self {
        Self {
            indices,
            rest_area,
            compliance,
            lambda: 0.0,
        }
    }
    /// Build an area constraint from the current configuration.
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
/// For two triangles sharing edge (p0,p1) with opposite vertices p2,p3 the
/// constraint value is `C = stiffness * sum_ab Q_ab p_a . p_b` with
/// `Q = 3 / (A0 + A1) * K K^T`, i.e. `C = 3 stiffness / (A0 + A1) * |sum_a K_a p_a|^2`.
/// `A0` and `A1` are the triangle areas and `K` holds cotangent weights that sum to
/// zero and annihilate the positions of every planar hinge, so `C` vanishes for a
/// flat hinge and is invariant under rigid translations and rotations.
#[derive(Debug, Clone)]
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
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Current constraint value `C = stiffness * sum_ab Q_ab p_a . p_b`.
    ///
    /// `C` is non-negative.  It is zero for a planar hinge and does not change
    /// when the four particles are translated or rotated rigidly.
    pub fn constraint_value(&self, particles: &[SoftParticle]) -> Real {
        let mut value = 0.0_f64;
        for (a, &ia) in self.indices.iter().enumerate() {
            for (b, &ib) in self.indices.iter().enumerate() {
                value +=
                    self.q_matrix[a * 4 + b] * particles[ia].position.dot(&particles[ib].position);
            }
        }
        value * self.stiffness
    }
    /// Cotangent weights `K` and the prefactor `3 / (A0 + A1)` of the hinge
    /// (p0, p1 | p2, p3), or `None` when both triangles are degenerate.
    ///
    /// With `c_p0(t)` and `c_p1(t)` the cotangents of the angles of triangle `t` at
    /// the shared-edge vertices p0 and p1 (`t1 = (p0, p1, p2)`, `t2 = (p0, p1, p3)`):
    ///
    /// `K = [c_p1(t1) + c_p1(t2), c_p0(t1) + c_p0(t2), -(c_p0(t1) + c_p1(t1)), -(c_p0(t2) + c_p1(t2))]`
    ///
    /// The angles at the wing vertices p2 and p3 do not enter.  `sum(K) = 0` makes
    /// `K . p` translation invariant, and `K . p = 0` for every planar hinge, so a
    /// flat sheet carries no bending energy.
    pub(super) fn hinge_weights(
        p0: &Vec3,
        p1: &Vec3,
        p2: &Vec3,
        p3: &Vec3,
    ) -> Option<([Real; 4], Real)> {
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
        // Angle at p0 / p1 in triangle (p0, p1, p2) and in triangle (p0, p1, p3).
        let cot_p0_t1 = cot(&e, &e0_2);
        let cot_p1_t1 = cot(&(-e), &e1_2);
        let cot_p0_t2 = cot(&e, &e0_3);
        let cot_p1_t2 = cot(&(-e), &e1_3);
        let a1 = 0.5 * e0_2.cross(&e).norm();
        let a2 = 0.5 * e0_3.cross(&e).norm();
        let a_total = a1 + a2;
        if a_total < 1e-14 {
            return None;
        }
        let k = [
            cot_p1_t1 + cot_p1_t2,
            cot_p0_t1 + cot_p0_t2,
            -(cot_p0_t1 + cot_p1_t1),
            -(cot_p0_t2 + cot_p1_t2),
        ];
        Some((k, 3.0 / a_total))
    }
    /// Quadratic form `Q = 3 / (A0 + A1) * K K^T` of the hinge (p0, p1 | p2, p3)
    /// (see [`hinge_weights`](Self::hinge_weights)); all zeros for a degenerate hinge.
    fn compute_q_matrix(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> [Real; 16] {
        let Some((k, scale)) = Self::hinge_weights(p0, p1, p2, p3) else {
            return [0.0; 16];
        };
        let mut q = [0.0_f64; 16];
        for row in 0..4 {
            for col in 0..4 {
                q[row * 4 + col] = scale * k[row] * k[col];
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
    /// Bergou quadratic form of the hinge; the same `Q` as
    /// [`IsometricBendingConstraint`] uses, so a flat hinge has zero energy.
    fn build_q(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> [Real; 16] {
        IsometricBendingConstraint::compute_q_matrix(p0, p1, p2, p3)
    }
}
/// Dihedral-angle bending constraint for two triangles sharing an edge.
///
/// Vertices (p0, p1) are the shared edge; p2 and p3 are the opposing vertices.
///
/// The constrained quantity is the angle `theta` in `[0, pi]` between the face
/// normals `(p2 - p0) x e` and `(p3 - p0) x e` with `e = p1 - p0`.  A flat hinge,
/// whose wings lie on opposite sides of the edge, has `theta = pi`; folding the
/// hinge in either direction lowers `theta`, and `theta = 0` means the two
/// triangles lie on top of each other.  The projection moves all four particles
/// along the gradient of `theta`, so it restores a fold in either direction and
/// conserves linear and angular momentum.
#[derive(Debug, Clone)]
pub struct BendingConstraint {
    /// Indices of the four particles: shared edge (i0, i1), wing vertices (i2, i3).
    pub indices: [usize; 4],
    /// Rest dihedral angle `theta` in radians, in `[0, pi]` (`pi` for a flat hinge).
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
    /// Compute the dihedral angle `theta` in `[0, pi]` between the two triangles
    /// sharing edge (p0,p1): `pi` for a flat hinge, smaller when it is folded.
    ///
    /// It is evaluated as `pi - |phi|` from the signed fold angle `phi`
    /// (`atan2`-based), which stays exact near the flat state where an `acos`
    /// of the normals' dot product would lose half of its digits.  A degenerate
    /// hinge (zero-length edge, wing on the hinge line) yields `0`.
    pub(super) fn compute_dihedral(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> Real {
        Self::hinge_fold(p0, p1, p2, p3).map_or(0.0, |fold| std::f64::consts::PI - fold.angle.abs())
    }
}
