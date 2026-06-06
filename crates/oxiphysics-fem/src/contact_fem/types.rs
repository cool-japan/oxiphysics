//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Mortar contact element for FEM — uses L2 projection of contact tractions.
///
/// The mortar method avoids locking and provides smooth stress fields.
pub struct MortarContactElement {
    /// Number of slave-side quadrature points.
    pub n_quad_pts: usize,
    /// Penalty parameter (or Lagrange multiplier stiffness).
    pub penalty: f64,
}
impl MortarContactElement {
    /// Compute mortar integral (L2-projected gap) over a 1D segment.
    ///
    /// Integrates the gap function weighted by shape functions using
    /// Gaussian quadrature (simplified 1D version).
    pub fn mortar_gap_integral(
        &self,
        slave_a: [f64; 3],
        slave_b: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
        surface_normal: [f64; 3],
    ) -> f64 {
        let (xi_pts, weights) = gauss_quad_1d(self.n_quad_pts);
        let mut integral = 0.0;
        let seg_len = {
            let d = [
                slave_b[0] - slave_a[0],
                slave_b[1] - slave_a[1],
                slave_b[2] - slave_a[2],
            ];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        };
        for (&xi, &w) in xi_pts.iter().zip(weights.iter()) {
            let t = 0.5 * (xi + 1.0);
            let slave_pt = [
                slave_a[0] + t * (slave_b[0] - slave_a[0]),
                slave_a[1] + t * (slave_b[1] - slave_a[1]),
                slave_a[2] + t * (slave_b[2] - slave_a[2]),
            ];
            let proj = project_point_to_segment(slave_pt, master_a, master_b);
            let d = [
                slave_pt[0] - proj[0],
                slave_pt[1] - proj[1],
                slave_pt[2] - proj[2],
            ];
            let gap =
                d[0] * surface_normal[0] + d[1] * surface_normal[1] + d[2] * surface_normal[2];
            integral += w * gap * seg_len * 0.5;
        }
        integral
    }
    /// Mortar constraint residual for a contact pair.
    ///
    /// R = penalty * integral(gap * N_I dS) for each slave node I.
    pub fn constraint_residual(
        &self,
        slave_a: [f64; 3],
        slave_b: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
        surface_normal: [f64; 3],
    ) -> f64 {
        let g_int = self.mortar_gap_integral(slave_a, slave_b, master_a, master_b, surface_normal);
        self.penalty * (-g_int).max(0.0)
    }
}
/// Stores candidate contact node pairs found by a proximity search.
pub struct ContactSearch {
    /// Pairs of node indices `(i_a, i_b)` that are within the search cutoff.
    pub candidate_pairs: Vec<(usize, usize)>,
}
impl ContactSearch {
    /// Find all node pairs from sets A and B whose distance is within `cutoff`.
    ///
    /// Complexity is O(|A| x |B|); suitable for small to medium meshes.
    pub fn brute_force_search(nodes_a: &[[f64; 3]], nodes_b: &[[f64; 3]], cutoff: f64) -> Self {
        let mut pairs = Vec::new();
        for (i, a) in nodes_a.iter().enumerate() {
            for (j, b) in nodes_b.iter().enumerate() {
                let dx = a[0] - b[0];
                let dy = a[1] - b[1];
                let dz = a[2] - b[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist <= cutoff {
                    pairs.push((i, j));
                }
            }
        }
        Self {
            candidate_pairs: pairs,
        }
    }
    /// Return the number of candidate contact pairs.
    pub fn pair_count(&self) -> usize {
        self.candidate_pairs.len()
    }
}
/// Penalty method for contact enforcement.
pub struct PenaltyContact {
    /// Penalty stiffness `k` (N/m).
    pub stiffness: f64,
    /// Coulomb friction coefficient `mu`.
    pub friction_coeff: f64,
}
impl PenaltyContact {
    /// Create a new [`PenaltyContact`].
    pub fn new(stiffness: f64, friction: f64) -> Self {
        Self {
            stiffness,
            friction_coeff: friction,
        }
    }
    /// Compute the normal contact force for the given gap.
    pub fn normal_force(&self, gap: f64) -> f64 {
        self.stiffness * (-gap).max(0.0)
    }
    /// Compute the tangential friction force using regularized Coulomb law.
    pub fn friction_force(&self, f_n: f64, sliding_vel: f64) -> f64 {
        const EPSILON: f64 = 1e-6;
        self.friction_coeff * f_n * (sliding_vel / EPSILON).tanh()
    }
    /// Compute both contact force components: `(normal_force, tangential_force)`.
    pub fn contact_residual(&self, gap: f64, sliding_vel: f64) -> (f64, f64) {
        let f_n = self.normal_force(gap);
        let f_t = self.friction_force(f_n, sliding_vel);
        (f_n, f_t)
    }
}
/// Augmented Lagrangian contact formulation.
pub struct AugmentedLagrangianContact {
    /// Penalty parameter `p`.
    pub penalty: f64,
    /// Normal Lagrange multiplier.
    pub lambda_n: f64,
    /// Tangential Lagrange multiplier.
    pub lambda_t: f64,
}
impl AugmentedLagrangianContact {
    /// Create a new [`AugmentedLagrangianContact`] with zero multipliers.
    pub fn new(penalty: f64) -> Self {
        Self {
            penalty,
            lambda_n: 0.0,
            lambda_t: 0.0,
        }
    }
    /// Update Lagrange multipliers from current gap and sliding displacement.
    pub fn update_multipliers(&mut self, gap: f64, sliding: f64, friction_coeff: f64) {
        self.lambda_n += self.penalty * gap.min(0.0);
        let trial = self.lambda_t + self.penalty * sliding;
        let limit = friction_coeff * self.lambda_n.abs();
        self.lambda_t = trial.clamp(-limit, limit);
    }
    /// Evaluate the augmented contact force for the current gap.
    pub fn force(&self, gap: f64) -> f64 {
        self.lambda_n + self.penalty * gap.min(0.0)
    }
}
/// Axis-Aligned Bounding Box for contact detection acceleration.
#[derive(Debug, Clone)]
pub struct Aabb {
    /// Minimum corner \[x, y, z\].
    pub min: [f64; 3],
    /// Maximum corner \[x, y, z\].
    pub max: [f64; 3],
}
impl Aabb {
    /// Create an AABB from a set of points.
    pub fn from_points(points: &[[f64; 3]]) -> Self {
        if points.is_empty() {
            return Self {
                min: [0.0; 3],
                max: [0.0; 3],
            };
        }
        let mut min = points[0];
        let mut max = points[0];
        for p in points.iter().skip(1) {
            for i in 0..3 {
                if p[i] < min[i] {
                    min[i] = p[i];
                }
                if p[i] > max[i] {
                    max[i] = p[i];
                }
            }
        }
        Self { min, max }
    }
    /// Expand the AABB by a margin in all directions.
    pub fn expand(&mut self, margin: f64) {
        for i in 0..3 {
            self.min[i] -= margin;
            self.max[i] += margin;
        }
    }
    /// Check if two AABBs overlap.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        for i in 0..3 {
            if self.max[i] < other.min[i] || self.min[i] > other.max[i] {
                return false;
            }
        }
        true
    }
    /// Check if a point is inside this AABB.
    pub fn contains_point(&self, point: [f64; 3]) -> bool {
        for (i, &pt) in point.iter().enumerate() {
            if pt < self.min[i] || pt > self.max[i] {
                return false;
            }
        }
        true
    }
    /// Compute the volume of the AABB.
    pub fn volume(&self) -> f64 {
        let mut vol = 1.0;
        for i in 0..3 {
            vol *= (self.max[i] - self.min[i]).max(0.0);
        }
        vol
    }
}
/// Accelerated contact search using AABB-based broad phase.
pub struct AcceleratedContactSearch;
impl AcceleratedContactSearch {
    /// Broad-phase contact detection using AABB overlap test.
    ///
    /// Groups nodes into element-level AABBs and only tests pairs
    /// whose AABBs overlap.
    ///
    /// Returns indices of element pairs that may be in contact.
    pub fn aabb_broad_phase(aabbs_a: &[Aabb], aabbs_b: &[Aabb]) -> Vec<(usize, usize)> {
        let mut candidate_pairs = Vec::new();
        for (i, a) in aabbs_a.iter().enumerate() {
            for (j, b) in aabbs_b.iter().enumerate() {
                if a.overlaps(b) {
                    candidate_pairs.push((i, j));
                }
            }
        }
        candidate_pairs
    }
    /// Bucket-based spatial hashing for contact detection.
    ///
    /// Assigns nodes to spatial buckets and only tests node pairs
    /// within the same or adjacent buckets.
    pub fn bucket_search(
        nodes_a: &[[f64; 3]],
        nodes_b: &[[f64; 3]],
        bucket_size: f64,
        cutoff: f64,
    ) -> Vec<(usize, usize)> {
        use std::collections::HashMap;
        let mut buckets: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for (j, b) in nodes_b.iter().enumerate() {
            let ix = (b[0] / bucket_size).floor() as i64;
            let iy = (b[1] / bucket_size).floor() as i64;
            let iz = (b[2] / bucket_size).floor() as i64;
            buckets.entry((ix, iy, iz)).or_default().push(j);
        }
        let mut pairs = Vec::new();
        for (i, a) in nodes_a.iter().enumerate() {
            let ix = (a[0] / bucket_size).floor() as i64;
            let iy = (a[1] / bucket_size).floor() as i64;
            let iz = (a[2] / bucket_size).floor() as i64;
            for dix in -1..=1 {
                for diy in -1..=1 {
                    for diz in -1..=1 {
                        let key = (ix + dix, iy + diy, iz + diz);
                        if let Some(indices) = buckets.get(&key) {
                            for &j in indices {
                                let b = &nodes_b[j];
                                let dx = a[0] - b[0];
                                let dy = a[1] - b[1];
                                let dz = a[2] - b[2];
                                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                                if dist <= cutoff {
                                    pairs.push((i, j));
                                }
                            }
                        }
                    }
                }
            }
        }
        pairs
    }
}
/// A node-to-segment contact element for 3D FEM.
///
/// Maps a slave node onto a master surface described by two end nodes,
/// computes gap function, normal, and penalty forces.
pub struct NodeToSegmentElement {
    /// Penalty stiffness for normal contact.
    pub penalty_normal: f64,
    /// Penalty stiffness for tangential contact.
    pub penalty_tangential: f64,
    /// Coulomb friction coefficient.
    pub friction_coeff: f64,
}
impl NodeToSegmentElement {
    /// Create a new node-to-segment element.
    pub fn new(penalty_normal: f64, penalty_tangential: f64, friction_coeff: f64) -> Self {
        Self {
            penalty_normal,
            penalty_tangential,
            friction_coeff,
        }
    }
    /// Compute signed gap from slave node to the nearest point on master segment.
    ///
    /// Positive gap → slave is outside (no contact).
    /// Negative gap → penetration.
    pub fn gap(
        slave_pos: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
        surface_normal: [f64; 3],
    ) -> f64 {
        let proj = project_point_to_segment(slave_pos, master_a, master_b);
        let d = [
            slave_pos[0] - proj[0],
            slave_pos[1] - proj[1],
            slave_pos[2] - proj[2],
        ];
        d[0] * surface_normal[0] + d[1] * surface_normal[1] + d[2] * surface_normal[2]
    }
    /// Compute nodal contact force vector for the slave node (penalty method).
    pub fn penalty_force(
        &self,
        slave_pos: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
        surface_normal: [f64; 3],
    ) -> [f64; 3] {
        let g = Self::gap(slave_pos, master_a, master_b, surface_normal);
        let f_n = self.penalty_normal * (-g).max(0.0);
        [
            f_n * surface_normal[0],
            f_n * surface_normal[1],
            f_n * surface_normal[2],
        ]
    }
    /// Shape function interpolation weight along segment at closest point.
    ///
    /// For segment \[a, b\] and closest point at parameter t ∈ \[0,1\]:
    ///   N_a = 1 - t, N_b = t.
    pub fn shape_functions(
        slave_pos: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
    ) -> [f64; 2] {
        let ab = [
            master_b[0] - master_a[0],
            master_b[1] - master_a[1],
            master_b[2] - master_a[2],
        ];
        let ap = [
            slave_pos[0] - master_a[0],
            slave_pos[1] - master_a[1],
            slave_pos[2] - master_a[2],
        ];
        let ab_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        if ab_sq < 1e-30 {
            return [1.0, 0.0];
        }
        let dot = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2];
        let t = (dot / ab_sq).clamp(0.0, 1.0);
        [1.0 - t, t]
    }
    /// Compute the 3x3 penalty contact stiffness matrix for slave node DOFs.
    pub fn contact_stiffness_matrix(
        &self,
        surface_normal: [f64; 3],
        in_contact: bool,
    ) -> [[f64; 3]; 3] {
        if !in_contact {
            return [[0.0; 3]; 3];
        }
        let k = self.penalty_normal;
        let mut ke = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                ke[i][j] = k * surface_normal[i] * surface_normal[j];
            }
        }
        ke
    }
}
/// Recover contact stress (pressure and shear) from nodal gap and forces.
pub struct ContactStressRecovery {
    /// Element area (m²) for stress recovery.
    pub element_area: f64,
}
impl ContactStressRecovery {
    /// Normal contact pressure from nodal force and element area.
    pub fn normal_pressure(&self, normal_force: f64) -> f64 {
        if self.element_area < 1e-30 {
            return 0.0;
        }
        normal_force / self.element_area
    }
    /// Shear stress from tangential force and element area.
    pub fn shear_stress(&self, tangential_force: f64) -> f64 {
        if self.element_area < 1e-30 {
            return 0.0;
        }
        tangential_force / self.element_area
    }
    /// Von Mises equivalent contact stress.
    pub fn von_mises_contact_stress(&self, normal_force: f64, tangential_force: f64) -> f64 {
        let sigma_n = self.normal_pressure(normal_force);
        let tau = self.shear_stress(tangential_force);
        (sigma_n * sigma_n + 3.0 * tau * tau).sqrt()
    }
    /// Maximum shear stress on contact interface.
    pub fn max_shear_stress(&self, normal_force: f64, tangential_force: f64) -> f64 {
        let sigma_n = self.normal_pressure(normal_force);
        let tau = self.shear_stress(tangential_force);
        ((sigma_n / 2.0).powi(2) + tau * tau).sqrt()
    }
    /// Compute nodal stress vector from element-level contact forces.
    ///
    /// Returns \[sigma_x, sigma_y, sigma_z\] using normal direction weighting.
    pub fn nodal_stress_vector(&self, normal_force: f64, normal_dir: [f64; 3]) -> [f64; 3] {
        let p = self.normal_pressure(normal_force);
        [p * normal_dir[0], p * normal_dir[1], p * normal_dir[2]]
    }
}
/// Tied contact constraint: enforces displacement continuity across an interface.
///
/// For tied contact, slave nodes are constrained to follow the master surface
/// (no relative motion allowed).
pub struct TiedContact {
    /// Pairs of (slave_node_index, master_node_index) that are tied.
    pub tied_pairs: Vec<(usize, usize)>,
    /// Penalty stiffness for enforcement.
    pub penalty: f64,
}
impl TiedContact {
    /// Create a new tied contact from node pairs.
    pub fn new(pairs: Vec<(usize, usize)>, penalty: f64) -> Self {
        Self {
            tied_pairs: pairs,
            penalty,
        }
    }
    /// Compute the tied contact force for a given pair.
    ///
    /// Force = penalty * (u_slave - u_master) applied to both nodes.
    pub fn tied_force(&self, slave_disp: [f64; 3], master_disp: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let mut f_slave = [0.0; 3];
        let mut f_master = [0.0; 3];
        for i in 0..3 {
            let gap = slave_disp[i] - master_disp[i];
            f_slave[i] = -self.penalty * gap;
            f_master[i] = self.penalty * gap;
        }
        (f_slave, f_master)
    }
    /// Compute the tied contact stiffness contribution for a pair.
    ///
    /// Returns the 6x6 stiffness matrix contribution (3 slave DOFs + 3 master DOFs).
    pub fn tied_stiffness(&self) -> [[f64; 6]; 6] {
        let k = self.penalty;
        let mut ke = [[0.0; 6]; 6];
        for i in 0..3 {
            ke[i][i] = k;
            ke[i][i + 3] = -k;
            ke[i + 3][i] = -k;
            ke[i + 3][i + 3] = k;
        }
        ke
    }
}
/// Sliding contact with Coulomb friction and stick-slip transition.
pub struct SlidingContact {
    /// Normal penalty stiffness.
    pub normal_stiffness: f64,
    /// Tangential penalty stiffness.
    pub tangential_stiffness: f64,
    /// Coulomb friction coefficient.
    pub friction_coeff: f64,
    /// Accumulated tangential displacement for stick detection.
    pub accumulated_slip: Vec<[f64; 3]>,
}
impl SlidingContact {
    /// Create a new sliding contact formulation.
    pub fn new(k_n: f64, k_t: f64, mu: f64, n_pairs: usize) -> Self {
        Self {
            normal_stiffness: k_n,
            tangential_stiffness: k_t,
            friction_coeff: mu,
            accumulated_slip: vec![[0.0; 3]; n_pairs],
        }
    }
    /// Compute the normal and tangential contact forces for a pair.
    ///
    /// Returns (force_normal_vec, force_tangential_vec, is_sliding).
    pub fn contact_forces(
        &self,
        gap: f64,
        normal: [f64; 3],
        tangential_disp: [f64; 3],
    ) -> ([f64; 3], [f64; 3], bool) {
        let pen = (-gap).max(0.0);
        let f_n = self.normal_stiffness * pen;
        let f_n_vec = [f_n * normal[0], f_n * normal[1], f_n * normal[2]];
        let mut f_t_trial = [0.0; 3];
        let mut f_t_mag_sq = 0.0;
        for i in 0..3 {
            f_t_trial[i] = self.tangential_stiffness * tangential_disp[i];
            f_t_mag_sq += f_t_trial[i] * f_t_trial[i];
        }
        let f_t_mag = f_t_mag_sq.sqrt();
        let f_t_limit = self.friction_coeff * f_n;
        let (f_t_vec, is_sliding);
        if f_t_mag <= f_t_limit || f_t_mag < 1e-30 {
            f_t_vec = f_t_trial;
            is_sliding = false;
        } else {
            let scale = f_t_limit / f_t_mag;
            f_t_vec = [
                f_t_trial[0] * scale,
                f_t_trial[1] * scale,
                f_t_trial[2] * scale,
            ];
            is_sliding = true;
        }
        (f_n_vec, f_t_vec, is_sliding)
    }
    /// Update accumulated slip for a contact pair.
    pub fn update_slip(&mut self, pair_idx: usize, slip_increment: [f64; 3]) {
        for (i, &inc) in slip_increment.iter().enumerate() {
            self.accumulated_slip[pair_idx][i] += inc;
        }
    }
    /// Reset accumulated slip for a pair (e.g., when contact is lost).
    pub fn reset_slip(&mut self, pair_idx: usize) {
        self.accumulated_slip[pair_idx] = [0.0; 3];
    }
}
