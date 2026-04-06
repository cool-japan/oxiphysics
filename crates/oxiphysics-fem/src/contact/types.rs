//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::f64::consts::PI;

/// A contact pair between two nodes with gap, normal, and multiplier.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactPair {
    /// Index of node A (slave).
    pub node_a: usize,
    /// Index of node B / face representative (master).
    pub node_b: usize,
    /// Gap function value (negative = penetration).
    pub gap: f64,
    /// Outward contact normal (pointing from master to slave).
    pub normal: [f64; 3],
    /// Lagrange multiplier.
    pub lambda: f64,
}
/// 2-D axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct Aabb2d {
    /// Minimum corner.
    pub min: [f64; 2],
    /// Maximum corner.
    pub max: [f64; 2],
}
impl Aabb2d {
    /// Create from min/max corners.
    #[allow(dead_code)]
    pub fn new(min: [f64; 2], max: [f64; 2]) -> Self {
        Self { min, max }
    }
    /// Expand box by `margin` on all sides.
    #[allow(dead_code)]
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            min: [self.min[0] - margin, self.min[1] - margin],
            max: [self.max[0] + margin, self.max[1] + margin],
        }
    }
    /// Test overlap with another AABB.
    #[allow(dead_code)]
    pub fn overlaps(&self, other: &Aabb2d) -> bool {
        self.max[0] >= other.min[0]
            && self.min[0] <= other.max[0]
            && self.max[1] >= other.min[1]
            && self.min[1] <= other.max[1]
    }
    /// Build from a list of 2-D points.
    #[allow(dead_code)]
    pub fn from_points(pts: &[[f64; 2]]) -> Option<Self> {
        if pts.is_empty() {
            return None;
        }
        let mut mn = pts[0];
        let mut mx = pts[0];
        for p in pts.iter().skip(1) {
            if p[0] < mn[0] {
                mn[0] = p[0];
            }
            if p[1] < mn[1] {
                mn[1] = p[1];
            }
            if p[0] > mx[0] {
                mx[0] = p[0];
            }
            if p[1] > mx[1] {
                mx[1] = p[1];
            }
        }
        Some(Self { min: mn, max: mx })
    }
    /// Center of the box.
    #[allow(dead_code)]
    pub fn center(&self) -> [f64; 2] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
        ]
    }
}
/// Node-to-segment (NTS) contact formulation.
///
/// One body is represented as discrete slave nodes; the other as master
/// segments.  The gap is computed by projecting each slave node onto the
/// closest master segment.
#[allow(dead_code)]
pub struct NodeToSegmentContact;
#[allow(dead_code)]
impl NodeToSegmentContact {
    /// Project a 3D slave node onto a master line segment and compute
    /// the signed gap.
    ///
    /// Returns `(xi, gap_vector)` where `xi ∈ [0,1]` is the parameter along
    /// the master segment and `gap_vector` is the shortest vector from the
    /// master segment to the slave node (positive = separated).
    pub fn project_node_to_segment_3d(
        slave: [f64; 3],
        seg_start: [f64; 3],
        seg_end: [f64; 3],
    ) -> (f64, [f64; 3]) {
        let ab = [
            seg_end[0] - seg_start[0],
            seg_end[1] - seg_start[1],
            seg_end[2] - seg_start[2],
        ];
        let ap = [
            slave[0] - seg_start[0],
            slave[1] - seg_start[1],
            slave[2] - seg_start[2],
        ];
        let ab_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        let xi = if ab_sq > 1e-30 {
            (ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / ab_sq
        } else {
            0.0
        }
        .clamp(0.0, 1.0);
        let closest = [
            seg_start[0] + xi * ab[0],
            seg_start[1] + xi * ab[1],
            seg_start[2] + xi * ab[2],
        ];
        let gap = [
            slave[0] - closest[0],
            slave[1] - closest[1],
            slave[2] - closest[2],
        ];
        (xi, gap)
    }
    /// Compute the signed gap (penetration negative) along the outward normal.
    pub fn signed_gap(
        slave: [f64; 3],
        seg_start: [f64; 3],
        seg_end: [f64; 3],
        outward_normal: [f64; 3],
    ) -> f64 {
        let (_, gap_vec) = Self::project_node_to_segment_3d(slave, seg_start, seg_end);
        gap_vec[0] * outward_normal[0]
            + gap_vec[1] * outward_normal[1]
            + gap_vec[2] * outward_normal[2]
    }
    /// Compute the contact force on the slave node using a penalty method.
    ///
    /// Returns the contact force vector (zero if no penetration).
    pub fn penalty_force(
        slave: [f64; 3],
        seg_start: [f64; 3],
        seg_end: [f64; 3],
        outward_normal: [f64; 3],
        penalty: f64,
    ) -> [f64; 3] {
        let gap = Self::signed_gap(slave, seg_start, seg_end, outward_normal);
        if gap >= 0.0 {
            return [0.0; 3];
        }
        let pen = -gap;
        [
            penalty * pen * outward_normal[0],
            penalty * pen * outward_normal[1],
            penalty * pen * outward_normal[2],
        ]
    }
}
/// Result of a Hertz contact calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct HertzResult {
    /// Contact radius `a` (m).
    pub contact_radius: f64,
    /// Normal contact force `F` (N).
    pub force: f64,
    /// Peak contact pressure `p0` (Pa).
    pub peak_pressure: f64,
    /// Approach (indentation depth) `delta` (m).
    pub approach: f64,
}
/// Utilities for computing the penetration depth (approach) from an applied force.
pub struct ContactPenetrationDepth;
impl ContactPenetrationDepth {
    /// Compute the indentation depth from an applied force for two elastic spheres.
    #[allow(clippy::too_many_arguments)]
    pub fn from_force_sphere_sphere(
        force: f64,
        r1: f64,
        r2: f64,
        e1: f64,
        nu1: f64,
        e2: f64,
        nu2: f64,
    ) -> f64 {
        let e_star = {
            let inv = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
            1.0 / inv
        };
        let r_star = 1.0 / (1.0 / r1 + 1.0 / r2);
        let numerator = 3.0 * force;
        let denominator = 4.0 * e_star * r_star.sqrt();
        (numerator / denominator).powf(2.0 / 3.0)
    }
}
/// Simple penalty contact enforcing non-penetration with Coulomb friction.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PenaltyContact {
    /// Normal penalty stiffness.
    pub penalty: f64,
    /// Coulomb friction coefficient.
    pub friction_coeff: f64,
}
#[allow(dead_code)]
impl PenaltyContact {
    /// Normal contact force: f = penalty * max(0, -gap).
    ///
    /// Positive when gap < 0 (penetration).
    pub fn normal_force(&self, gap: f64) -> f64 {
        self.penalty * (-gap).max(0.0)
    }
    /// Tangential (friction) force via Coulomb law.
    ///
    /// `slip` is the tangential slip vector; `normal_f` is the magnitude of
    /// the normal contact force.  Returns the friction force vector opposing
    /// slip, bounded by μ * |normal_f|.
    pub fn tangential_force(&self, slip: [f64; 3], normal_f: f64) -> [f64; 3] {
        let slip_norm = (slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2]).sqrt();
        if slip_norm < 1e-30 || normal_f <= 0.0 {
            return [0.0; 3];
        }
        let limit = self.friction_coeff * normal_f;
        let scale = limit / slip_norm;
        [-slip[0] * scale, -slip[1] * scale, -slip[2] * scale]
    }
}
/// Mortar-based contact formulation.
///
/// Uses weighted integrals of shape functions on the contact surface to
/// enforce contact constraints in a variationally consistent manner.
#[allow(dead_code)]
pub struct MortarContact;
impl MortarContact {
    /// Compute the mortar integral (D matrix entry) between a slave node
    /// and a master segment.
    ///
    /// D_ij = integral over master segment of N_slave_i * N_master_j ds
    ///
    /// Simplified: uses trapezoidal rule on a 1D segment.
    #[allow(dead_code)]
    pub fn mortar_d_integral(slave_xi: f64, segment_length: f64) -> [f64; 2] {
        let n1 = 1.0 - slave_xi;
        let n2 = slave_xi;
        [n1 * segment_length, n2 * segment_length]
    }
    /// Compute the mortar mass matrix M entry.
    ///
    /// M_ij = integral over slave segment of N_slave_i * N_slave_j ds
    #[allow(dead_code)]
    pub fn mortar_m_integral(segment_length: f64) -> [[f64; 2]; 2] {
        let l = segment_length;
        [[2.0 * l / 6.0, l / 6.0], [l / 6.0, 2.0 * l / 6.0]]
    }
    /// Project a slave node onto a master surface (1D case).
    ///
    /// Given slave node position and master segment endpoints,
    /// returns the parametric coordinate xi in \[0,1\] on the master segment
    /// and the gap distance.
    #[allow(dead_code)]
    pub fn project_slave_to_master(
        slave_pos: [f64; 2],
        master_start: [f64; 2],
        master_end: [f64; 2],
    ) -> (f64, f64) {
        let dx = master_end[0] - master_start[0];
        let dy = master_end[1] - master_start[1];
        let seg_len_sq = dx * dx + dy * dy;
        if seg_len_sq < 1e-30 {
            let gap_x = slave_pos[0] - master_start[0];
            let gap_y = slave_pos[1] - master_start[1];
            return (0.0, (gap_x * gap_x + gap_y * gap_y).sqrt());
        }
        let t_x = slave_pos[0] - master_start[0];
        let t_y = slave_pos[1] - master_start[1];
        let xi = ((t_x * dx + t_y * dy) / seg_len_sq).clamp(0.0, 1.0);
        let cp = [master_start[0] + xi * dx, master_start[1] + xi * dy];
        let gap_x = slave_pos[0] - cp[0];
        let gap_y = slave_pos[1] - cp[1];
        let gap = (gap_x * gap_x + gap_y * gap_y).sqrt();
        (xi, gap)
    }
}
/// Axis-aligned bounding box (AABB) for contact detection.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
#[allow(dead_code)]
impl Aabb {
    /// Create a new AABB.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Expand by a margin in all directions.
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }
    /// Check if this AABB overlaps with another.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Check if a point is inside this AABB.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        point[0] >= self.min[0]
            && point[0] <= self.max[0]
            && point[1] >= self.min[1]
            && point[1] <= self.max[1]
            && point[2] >= self.min[2]
            && point[2] <= self.max[2]
    }
    /// Compute AABB from a list of points.
    pub fn from_points(points: &[[f64; 3]]) -> Option<Self> {
        if points.is_empty() {
            return None;
        }
        let mut min = points[0];
        let mut max = points[0];
        for p in points.iter().skip(1) {
            for k in 0..3 {
                if p[k] < min[k] {
                    min[k] = p[k];
                }
                if p[k] > max[k] {
                    max[k] = p[k];
                }
            }
        }
        Some(Self { min, max })
    }
}
/// Assemble contact stiffness contributions as sparse (row, col, value) triplets.
///
/// For each contact pair the diagonal DOF contributions are added.
#[allow(dead_code)]
pub struct ContactStiffness;
#[allow(dead_code)]
impl ContactStiffness {
    /// Assemble contact stiffness triplets.
    ///
    /// Each active contact pair contributes penalty to the diagonal DOFs of
    /// `node_a` and `node_b` (3 DOF each, all three spatial directions).
    /// Returns a list of `(row, col, value)` triplets.
    pub fn assemble(
        pairs: &[ContactPair],
        penalty: f64,
        _n_dof: usize,
    ) -> Vec<(usize, usize, f64)> {
        let mut triplets = Vec::new();
        for pair in pairs {
            if pair.gap >= 0.0 {
                continue;
            }
            for k in 0..3 {
                let dof_a = 3 * pair.node_a + k;
                let dof_b = 3 * pair.node_b + k;
                let kc = penalty * pair.normal[k] * pair.normal[k];
                triplets.push((dof_a, dof_a, kc));
                triplets.push((dof_b, dof_b, kc));
                triplets.push((dof_a, dof_b, -kc));
                triplets.push((dof_b, dof_a, -kc));
            }
        }
        triplets
    }
}
/// Validate the Hertz contact theory predictions by checking the relationship
/// between approach, contact radius, and force.
///
/// For two identical spheres:
/// - a = sqrt(R* * δ)
/// - F = 4/3 * E* * R*^(1/2) * δ^(3/2)
/// - p0 = 3F / (2π a²) = (6 F E*² / (π³ R*²))^(1/3)
#[allow(dead_code)]
pub struct HertzValidator;
#[allow(dead_code)]
impl HertzValidator {
    /// Check the Hertz force-approach relationship: F = 4/3 E* sqrt(R) δ^1.5
    ///
    /// Returns relative error between computed and expected force.
    pub fn validate_force(result: &HertzResult, e_star: f64, r_star: f64) -> f64 {
        let expected = (4.0 / 3.0) * e_star * r_star.sqrt() * result.approach.powf(1.5);
        if expected.abs() < 1e-30 {
            return 0.0;
        }
        (result.force - expected).abs() / expected
    }
    /// Check the contact radius formula: a² = R* δ
    pub fn validate_contact_radius(result: &HertzResult, r_star: f64) -> f64 {
        let expected_a = (r_star * result.approach).sqrt();
        if expected_a.abs() < 1e-30 {
            return 0.0;
        }
        (result.contact_radius - expected_a).abs() / expected_a
    }
    /// Check the peak pressure: p0 = 3F/(2πa²)
    pub fn validate_peak_pressure(result: &HertzResult) -> f64 {
        if result.contact_radius < 1e-30 {
            return 0.0;
        }
        let expected = 3.0 * result.force / (2.0 * PI * result.contact_radius.powi(2));
        if expected.abs() < 1e-30 {
            return 0.0;
        }
        (result.peak_pressure - expected).abs() / expected
    }
    /// Compute the mean contact pressure: p_mean = F / (π a²)
    pub fn mean_pressure(result: &HertzResult) -> f64 {
        if result.contact_radius < 1e-30 {
            return 0.0;
        }
        result.force / (PI * result.contact_radius.powi(2))
    }
    /// Check that peak pressure = 1.5 * mean pressure (Hertz identity)
    pub fn validate_pressure_ratio(result: &HertzResult) -> f64 {
        let p_mean = Self::mean_pressure(result);
        if p_mean.abs() < 1e-30 {
            return 0.0;
        }
        let expected_ratio = 1.5;
        let actual_ratio = result.peak_pressure / p_mean;
        (actual_ratio - expected_ratio).abs() / expected_ratio
    }
}
/// Contact detection between two FEM node sets using AABB trees.
#[allow(dead_code)]
pub struct FemContactDetector {
    /// Skin distance: a node is a candidate if it is within this distance of the other body.
    pub skin_distance: f64,
}
#[allow(dead_code)]
impl FemContactDetector {
    /// Create a new contact detector.
    pub fn new(skin_distance: f64) -> Self {
        Self { skin_distance }
    }
    /// Find candidate contact pairs between two node sets.
    ///
    /// Returns a list of `(slave_node_idx, master_node_idx)` pairs
    /// whose distance is within `skin_distance`.
    pub fn find_candidates(
        &self,
        slave_nodes: &[[f64; 3]],
        master_nodes: &[[f64; 3]],
    ) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for (i, s) in slave_nodes.iter().enumerate() {
            for (j, m) in master_nodes.iter().enumerate() {
                let dx = s[0] - m[0];
                let dy = s[1] - m[1];
                let dz = s[2] - m[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist <= self.skin_distance {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
    /// Find all slave nodes within the master body AABB.
    pub fn nodes_inside_aabb(&self, slave_nodes: &[[f64; 3]], master_aabb: &Aabb) -> Vec<usize> {
        let expanded = master_aabb.expanded(self.skin_distance);
        slave_nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| expanded.contains(**p))
            .map(|(i, _)| i)
            .collect()
    }
}
/// Coulomb friction status of a contact pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum FrictionStatus {
    /// No contact (gap >= 0).
    Open,
    /// Contact, no sliding (slip < tol).
    Stick,
    /// Contact, sliding.
    Slip,
}
/// Elliptical Hertz contact solver (approximate — circular limit).
///
/// For a general elliptic contact problem, the contact patch is elliptic.
/// Here we use the simplification that for principal curvatures `kappa1` and `kappa2`
/// the effective radius is `R* = sqrt(R_x * R_y)` (geometric mean).
#[allow(dead_code)]
pub struct EllipticalHertz;
impl EllipticalHertz {
    /// Approximate elliptical contact from two equivalent radii.
    ///
    /// `R_x` and `R_y` are the principal equivalent radii.
    /// Result semi-axes are estimated assuming `a/b = (R_x/R_y)^(1/3)`.
    pub fn contact(r_x: f64, r_y: f64, e_star: f64, force: f64) -> EllipticalHertzResult {
        let r_eff = (r_x * r_y).sqrt();
        let delta = (force * force / (e_star * e_star * r_eff)).cbrt() * 0.5;
        let a_circ = (r_eff * delta).sqrt();
        let ratio = (r_x / r_y).powf(1.0 / 3.0);
        let a = a_circ * ratio.sqrt();
        let b = a_circ / ratio.sqrt();
        let p0 = 3.0 * force / (2.0 * PI * a * b);
        EllipticalHertzResult {
            semi_axis_a: a,
            semi_axis_b: b,
            peak_pressure: p0,
            approach: delta,
        }
    }
    /// Contact pressure distribution over elliptic patch.
    ///
    /// `p(x,y) = p0 * sqrt(1 - (x/a)^2 - (y/b)^2)` for (x/a)^2+(y/b)^2 <= 1, else 0.
    pub fn pressure_distribution(x: f64, y: f64, p0: f64, semi_a: f64, semi_b: f64) -> f64 {
        let arg = (x / semi_a) * (x / semi_a) + (y / semi_b) * (y / semi_b);
        if arg > 1.0 {
            0.0
        } else {
            p0 * (1.0 - arg).sqrt()
        }
    }
}
/// Hertz contact mechanics solver for common contact geometries.
pub struct HertzContact;
impl HertzContact {
    /// Compute reduced (combined) elastic modulus `E*`.
    #[inline]
    fn reduced_modulus(e1: f64, nu1: f64, e2: f64, nu2: f64) -> f64 {
        let inv = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
        1.0 / inv
    }
    /// Solve Hertz contact for two elastic spheres pressed together.
    #[allow(clippy::too_many_arguments)]
    pub fn sphere_sphere(
        r1: f64,
        r2: f64,
        e1: f64,
        nu1: f64,
        e2: f64,
        nu2: f64,
        approach: f64,
    ) -> HertzResult {
        let e_star = Self::reduced_modulus(e1, nu1, e2, nu2);
        let r_star = 1.0 / (1.0 / r1 + 1.0 / r2);
        let a = (r_star * approach).sqrt();
        let force = (4.0 / 3.0) * e_star * r_star.sqrt() * approach.powf(1.5);
        let peak_pressure = if a > 0.0 {
            3.0 * force / (2.0 * PI * a * a)
        } else {
            0.0
        };
        HertzResult {
            contact_radius: a,
            force,
            peak_pressure,
            approach,
        }
    }
    /// Solve Hertz contact for an elastic sphere on a flat surface.
    #[allow(clippy::too_many_arguments)]
    pub fn sphere_flat(r: f64, e1: f64, nu1: f64, e2: f64, nu2: f64, approach: f64) -> HertzResult {
        let e_star = Self::reduced_modulus(e1, nu1, e2, nu2);
        let r_star = r;
        let a = (r_star * approach).sqrt();
        let force = (4.0 / 3.0) * e_star * r_star.sqrt() * approach.powf(1.5);
        let peak_pressure = if a > 0.0 {
            3.0 * force / (2.0 * PI * a * a)
        } else {
            0.0
        };
        HertzResult {
            contact_radius: a,
            force,
            peak_pressure,
            approach,
        }
    }
    /// Solve 2D Hertz contact for a cylinder pressed against a flat surface.
    #[allow(clippy::too_many_arguments)]
    pub fn cylinder_flat_2d(
        r: f64,
        e1: f64,
        nu1: f64,
        e2: f64,
        nu2: f64,
        half_length: f64,
        approach: f64,
    ) -> HertzResult {
        let e_star = Self::reduced_modulus(e1, nu1, e2, nu2);
        let r_star = r;
        let b = (4.0 * r_star * approach / PI).sqrt();
        let length = 2.0 * half_length;
        let force = (PI / 2.0) * e_star * approach * length;
        let peak_pressure = if b > 0.0 && length > 0.0 {
            2.0 * force / (PI * b * length)
        } else {
            0.0
        };
        HertzResult {
            contact_radius: b,
            force,
            peak_pressure,
            approach,
        }
    }
}
/// Dual (multiplier) state for augmented Lagrangian frictionless contact.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DualContactState {
    /// Lagrange multipliers λ (contact pressure, one per pair).
    pub lambda: Vec<f64>,
    /// Penalty parameter `r`.
    pub r: f64,
}
impl DualContactState {
    /// Initialise with zero multipliers.
    #[allow(dead_code)]
    pub fn new(n_pairs: usize, r: f64) -> Self {
        Self {
            lambda: vec![0.0; n_pairs],
            r,
        }
    }
    /// Uzawa update step: `lambda_{k+1} = max(0, lambda_k + r * g_n)`.
    #[allow(dead_code)]
    pub fn uzawa_step(&mut self, gaps: &[f64]) {
        for (lam, &g) in self.lambda.iter_mut().zip(gaps.iter()) {
            *lam = (*lam + self.r * g).min(0.0);
        }
    }
    /// Active set: indices where contact is active (lambda < 0).
    #[allow(dead_code)]
    pub fn active_set(&self) -> Vec<usize> {
        self.lambda
            .iter()
            .enumerate()
            .filter(|&(_, l)| *l < 0.0)
            .map(|(i, _)| i)
            .collect()
    }
    /// Contact force contribution `f_c = -lambda * n` at each pair.
    #[allow(dead_code)]
    pub fn contact_forces(&self, normals: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.lambda
            .iter()
            .zip(normals.iter())
            .map(|(&lam, n)| [-lam * n[0], -lam * n[1], -lam * n[2]])
            .collect()
    }
}
/// Result of an elliptical Hertz contact calculation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct EllipticalHertzResult {
    /// Semi-axis in x-direction `a` (m).
    pub semi_axis_a: f64,
    /// Semi-axis in y-direction `b` (m).
    pub semi_axis_b: f64,
    /// Peak contact pressure `p0` (Pa).
    pub peak_pressure: f64,
    /// Approach (indentation) `delta` (m).
    pub approach: f64,
}
/// Augmented Lagrangian contact formulation with Uzawa update.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AugmentedLagrangianSolver {
    /// Penalty parameter.
    pub penalty: f64,
    /// Normal Lagrange multipliers (one per contact pair).
    pub lambda_n: Vec<f64>,
    /// Tangential Lagrange multipliers (one per contact pair).
    pub lambda_t: Vec<f64>,
    /// Uzawa update tolerance.
    pub tol: f64,
    /// Maximum Uzawa iterations.
    pub max_iter: usize,
}
impl AugmentedLagrangianSolver {
    /// Create a new augmented Lagrangian solver for n_pairs contact pairs.
    #[allow(dead_code)]
    pub fn new(penalty: f64, n_pairs: usize, tol: f64, max_iter: usize) -> Self {
        Self {
            penalty,
            lambda_n: vec![0.0; n_pairs],
            lambda_t: vec![0.0; n_pairs],
            tol,
            max_iter,
        }
    }
    /// Perform Uzawa update for all contact pairs.
    ///
    /// Returns true if converged (max multiplier change < tol).
    #[allow(dead_code)]
    pub fn uzawa_update(&mut self, gaps: &[f64], slidings: &[f64], friction_coeff: f64) -> bool {
        assert_eq!(gaps.len(), self.lambda_n.len());
        assert_eq!(slidings.len(), self.lambda_t.len());
        let mut max_change = 0.0_f64;
        for i in 0..self.lambda_n.len() {
            let old_n = self.lambda_n[i];
            self.lambda_n[i] = (self.lambda_n[i] + self.penalty * gaps[i]).min(0.0);
            let old_t = self.lambda_t[i];
            let trial = self.lambda_t[i] + self.penalty * slidings[i];
            let limit = friction_coeff * self.lambda_n[i].abs();
            self.lambda_t[i] = trial.clamp(-limit, limit);
            max_change = max_change.max((self.lambda_n[i] - old_n).abs());
            max_change = max_change.max((self.lambda_t[i] - old_t).abs());
        }
        max_change < self.tol
    }
    /// Compute the augmented normal force for a contact pair.
    #[allow(dead_code)]
    pub fn normal_force(&self, pair: usize, gap: f64) -> f64 {
        (self.lambda_n[pair] + self.penalty * gap).min(0.0)
    }
}
/// Configuration for penalty-based contact enforcement.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PenaltyParameters {
    /// Normal penalty stiffness (N/m).
    pub normal_stiffness: f64,
    /// Tangential penalty stiffness (N/m).
    pub tangential_stiffness: f64,
    /// Regularization velocity for friction (m/s).
    pub regularization_velocity: f64,
    /// Maximum allowed penetration before stiffness scaling (m).
    pub max_penetration: f64,
    /// Stiffness scaling factor for excessive penetration.
    pub scaling_factor: f64,
}
impl PenaltyParameters {
    /// Create default penalty parameters based on material stiffness.
    ///
    /// A common heuristic: penalty stiffness = alpha * E * A / h
    /// where alpha is a scaling factor (typically 10-1000).
    #[allow(dead_code)]
    pub fn from_material(youngs_modulus: f64, element_size: f64, alpha: f64) -> Self {
        let k_n = alpha * youngs_modulus / element_size;
        Self {
            normal_stiffness: k_n,
            tangential_stiffness: k_n * 0.5,
            regularization_velocity: 1e-6,
            max_penetration: element_size * 0.1,
            scaling_factor: 10.0,
        }
    }
    /// Compute the adaptive normal penalty stiffness based on current penetration.
    ///
    /// If penetration exceeds max_penetration, the stiffness is increased.
    #[allow(dead_code)]
    pub fn adaptive_stiffness(&self, penetration: f64) -> f64 {
        if penetration > self.max_penetration {
            self.normal_stiffness * self.scaling_factor
        } else {
            self.normal_stiffness
        }
    }
    /// Compute penalty contact force (normal component).
    #[allow(dead_code)]
    pub fn normal_force(&self, gap: f64) -> f64 {
        let pen = (-gap).max(0.0);
        let k = self.adaptive_stiffness(pen);
        k * pen
    }
    /// Compute regularized friction force.
    #[allow(dead_code)]
    pub fn friction_force(&self, normal_force: f64, sliding_vel: f64, friction_coeff: f64) -> f64 {
        friction_coeff * normal_force * (sliding_vel / self.regularization_velocity).tanh()
    }
}
/// Augmented Lagrangian contact with per-contact-point multipliers.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AugmentedLagrangianContact {
    /// Penalty parameter.
    pub penalty: f64,
    /// Lagrange multipliers (one per contact point).
    pub lambda: Vec<f64>,
}
#[allow(dead_code)]
impl AugmentedLagrangianContact {
    /// Create a new AL contact object with `n` zero multipliers.
    pub fn new(penalty: f64, n: usize) -> Self {
        Self {
            penalty,
            lambda: vec![0.0; n],
        }
    }
    /// Uzawa multiplier update: λ_i ← max(0, λ_i + penalty * gap_i).
    ///
    /// For contact, gap < 0 means penetration.  Active contact increases λ.
    pub fn update_multipliers(&mut self, gaps: &[f64]) {
        assert_eq!(gaps.len(), self.lambda.len());
        for (lam, &g) in self.lambda.iter_mut().zip(gaps.iter()) {
            *lam = (*lam - self.penalty * g).max(0.0);
        }
    }
}
/// Segment-to-segment contact formulation for 2D contact.
#[allow(dead_code)]
pub struct SegmentToSegmentContact;
impl SegmentToSegmentContact {
    /// Compute the minimum distance between two line segments in 2D.
    ///
    /// Segment 1: p1_start to p1_end
    /// Segment 2: p2_start to p2_end
    ///
    /// Returns (distance, parameter_t1, parameter_t2) where t1, t2 in \[0,1\]
    /// parameterize the closest points on each segment.
    #[allow(dead_code)]
    pub fn segment_distance(
        p1_start: [f64; 2],
        p1_end: [f64; 2],
        p2_start: [f64; 2],
        p2_end: [f64; 2],
    ) -> (f64, f64, f64) {
        let d1 = [p1_end[0] - p1_start[0], p1_end[1] - p1_start[1]];
        let d2 = [p2_end[0] - p2_start[0], p2_end[1] - p2_start[1]];
        let r = [p1_start[0] - p2_start[0], p1_start[1] - p2_start[1]];
        let a = d1[0] * d1[0] + d1[1] * d1[1];
        let e = d2[0] * d2[0] + d2[1] * d2[1];
        let f = d2[0] * r[0] + d2[1] * r[1];
        let b = d1[0] * d2[0] + d1[1] * d2[1];
        let c = d1[0] * r[0] + d1[1] * r[1];
        let denom = a * e - b * b;
        let (t1, t2);
        if denom.abs() < 1e-30 {
            t1 = 0.0;
            t2 = if e.abs() > 1e-30 { f / e } else { 0.0 };
        } else {
            t1 = ((b * f - c * e) / denom).clamp(0.0, 1.0);
            t2 = ((a * f - b * c) / denom).clamp(0.0, 1.0);
        }
        let cp1 = [p1_start[0] + t1 * d1[0], p1_start[1] + t1 * d1[1]];
        let cp2 = [p2_start[0] + t2 * d2[0], p2_start[1] + t2 * d2[1]];
        let diff = [cp1[0] - cp2[0], cp1[1] - cp2[1]];
        let dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();
        (dist, t1, t2)
    }
    /// Compute the gap function for segment-to-segment contact.
    ///
    /// Positive gap means separation; negative means penetration.
    #[allow(dead_code)]
    pub fn gap_function(
        p1_start: [f64; 2],
        p1_end: [f64; 2],
        p2_start: [f64; 2],
        p2_end: [f64; 2],
        outward_normal: [f64; 2],
    ) -> f64 {
        let (dist, t1, _t2) = Self::segment_distance(p1_start, p1_end, p2_start, p2_end);
        let d1 = [p1_end[0] - p1_start[0], p1_end[1] - p1_start[1]];
        let d2 = [p2_end[0] - p2_start[0], p2_end[1] - p2_start[1]];
        let cp1 = [p1_start[0] + t1 * d1[0], p1_start[1] + t1 * d1[1]];
        let _t2_val = _t2;
        let cp2 = [p2_start[0] + _t2_val * d2[0], p2_start[1] + _t2_val * d2[1]];
        let diff = [cp1[0] - cp2[0], cp1[1] - cp2[1]];
        let signed_gap = diff[0] * outward_normal[0] + diff[1] * outward_normal[1];
        if signed_gap.abs() < 1e-30 {
            dist
        } else {
            signed_gap
        }
    }
}
