//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use std::f64::consts::PI;

/// Coupling zone between a standard FEM region and an XFEM enriched region.
///
/// A ramp (blending) function is applied over elements in the coupling zone to
/// ensure a smooth transition of the enrichment.
#[derive(Debug, Clone)]
pub struct XfemFemCouplingZone {
    /// Node indices in the blending/coupling layer.
    pub blend_nodes: Vec<usize>,
    /// Ramp weights for each blending node (0 = pure FEM, 1 = pure XFEM).
    pub ramp_weights: Vec<f64>,
}
impl XfemFemCouplingZone {
    /// Create a coupling zone by linear ramping from 0 to 1 over `n_layers`
    /// of nodes starting at `start_node_id`.
    pub fn linear_ramp(start_node_id: usize, n_layers: usize) -> Self {
        let blend_nodes: Vec<usize> = (start_node_id..start_node_id + n_layers).collect();
        let ramp_weights: Vec<f64> = (0..n_layers)
            .map(|i| i as f64 / (n_layers - 1).max(1) as f64)
            .collect();
        Self {
            blend_nodes,
            ramp_weights,
        }
    }
    /// Evaluate the ramp function at a node by its position in the blending list.
    pub fn ramp_at(&self, local_idx: usize) -> f64 {
        self.ramp_weights.get(local_idx).copied().unwrap_or(1.0)
    }
    /// Blended enrichment value: `w * xfem_val + (1 - w) * fem_val`.
    pub fn blend(&self, local_idx: usize, xfem_val: f64, fem_val: f64) -> f64 {
        let w = self.ramp_at(local_idx);
        w * xfem_val + (1.0 - w) * fem_val
    }
}
/// Blending (transition/ramp) function for XFEM enrichment.
///
/// In elements that are partially enriched (some nodes enriched, some not),
/// a blending function ensures partition of unity.
///
/// The ramp function R(x) = sum_J N_J(x) where J ranges over enriched nodes.
pub struct BlendingFunction;
impl BlendingFunction {
    /// Compute the ramp function value for a partially enriched element.
    ///
    /// `shape_functions` - standard FEM shape functions N_I at the evaluation point
    /// `enriched_flags` - boolean flag for each node (true if enriched)
    ///
    /// Returns the ramp function R = sum over enriched nodes of N_I.
    pub fn ramp(shape_functions: &[f64], enriched_flags: &[bool]) -> f64 {
        assert_eq!(shape_functions.len(), enriched_flags.len());
        let mut r = 0.0;
        for (n, &is_enriched) in shape_functions.iter().zip(enriched_flags.iter()) {
            if is_enriched {
                r += n;
            }
        }
        r
    }
    /// Compute the modified enrichment function for a blending element.
    ///
    /// The enrichment function is multiplied by the ramp function to ensure
    /// partition of unity: F_modified = R(x) * F(x) for fully enriched,
    /// or F_modified = N_I * (F(x) - F(x_I)) for shifted enrichment.
    ///
    /// Returns the shifted enrichment values for each enriched node.
    pub fn shifted_enrichment(enrichment_values: &[f64], enrichment_at_nodes: &[f64]) -> Vec<f64> {
        assert_eq!(enrichment_values.len(), enrichment_at_nodes.len());
        enrichment_values
            .iter()
            .zip(enrichment_at_nodes.iter())
            .map(|(&f, &f_node)| f - f_node)
            .collect()
    }
}
/// An enriched node carrying extra degrees of freedom for XFEM.
pub struct EnrichedNode {
    /// Index of the underlying FEM node.
    pub node_id: usize,
    /// Which enrichment type is applied.
    pub enrichment: EnrichmentType,
    /// Extra enrichment DOF values.
    pub dofs: Vec<f64>,
}
/// Utilities for computing and using Stress Intensity Factors.
pub struct SifCalculator {
    /// Young's modulus \[Pa\].
    pub e_modulus: f64,
    /// Poisson's ratio \[dimensionless\].
    pub poisson: f64,
}
impl SifCalculator {
    /// Construct a new `SifCalculator`.
    pub fn new(e: f64, nu: f64) -> Self {
        SifCalculator {
            e_modulus: e,
            poisson: nu,
        }
    }
    /// Approximate mode-I stress intensity factor from crack-opening displacement.
    pub fn ki_from_displacement(crack: &Crack, u_y: f64, r: f64) -> f64 {
        let _ = crack;
        u_y * (2.0 * PI / r).sqrt()
    }
    /// Irwin plastic-zone size estimate.
    pub fn plastic_zone_size(&self, k1: f64, yield_stress: f64) -> f64 {
        let denom = yield_stress * (2.0 * PI).sqrt();
        (k1 / denom).powi(2)
    }
}
/// Level set representation for the crack geometry.
///
/// Uses two level set functions:
/// - phi: signed distance to the crack surface (zero on the crack plane)
/// - psi: signed distance to the crack front/tip (zero at the tip)
///
/// A point is on the crack if phi = 0 and psi < 0.
pub struct LevelSet {
    /// Nodal values of the normal level set phi (signed distance to crack plane).
    pub phi: Vec<f64>,
    /// Nodal values of the tangential level set psi (signed distance to crack front).
    pub psi: Vec<f64>,
}
impl LevelSet {
    /// Create level sets for a straight crack from `start` to `tip` for a given set of nodes.
    pub fn from_crack(nodes: &[[f64; 3]], crack: &Crack) -> Self {
        let last = crack.segments.last().expect("crack must have segments");
        let n = nodes.len();
        let mut phi = vec![0.0; n];
        let mut psi = vec![0.0; n];
        for (i, node) in nodes.iter().enumerate() {
            phi[i] = signed_distance_to_line(*node, last.start, last.end);
            let tip = crack.tip;
            let dir = crack.growth_direction;
            let dx = node[0] - tip[0];
            let dy = node[1] - tip[1];
            let dz = node[2] - tip[2];
            psi[i] = dx * dir[0] + dy * dir[1] + dz * dir[2];
        }
        LevelSet { phi, psi }
    }
    /// Determine if a node is enriched with Heaviside enrichment.
    ///
    /// A node needs Heaviside enrichment if the element it belongs to
    /// is cut by the crack (phi changes sign across the element).
    /// Simplified check: node has psi < 0 (behind the crack tip).
    pub fn needs_heaviside(&self, node_idx: usize) -> bool {
        self.psi[node_idx] < 0.0
    }
    /// Determine if a node is enriched with crack-tip enrichment.
    ///
    /// A node needs tip enrichment if it is near the crack tip.
    /// Simplified check: |psi| < enrichment_radius and |phi| < enrichment_radius.
    pub fn needs_tip_enrichment(&self, node_idx: usize, enrichment_radius: f64) -> bool {
        self.psi[node_idx].abs() < enrichment_radius && self.phi[node_idx].abs() < enrichment_radius
    }
    /// Update level sets after crack propagation by `da` in `direction`.
    pub fn update_after_propagation(
        &mut self,
        nodes: &[[f64; 3]],
        new_tip: [f64; 3],
        new_direction: [f64; 3],
    ) {
        for (i, node) in nodes.iter().enumerate() {
            let dx = node[0] - new_tip[0];
            let dy = node[1] - new_tip[1];
            let dz = node[2] - new_tip[2];
            self.psi[i] = dx * new_direction[0] + dy * new_direction[1] + dz * new_direction[2];
        }
    }
    /// Interpolate the level set value at an arbitrary point using
    /// linear interpolation between two nodes.
    pub fn interpolate_phi(&self, node_a: usize, node_b: usize, xi: f64) -> f64 {
        (1.0 - xi) * self.phi[node_a] + xi * self.phi[node_b]
    }
    /// Find the zero-crossing parameter xi in \[0,1\] between nodes a and b.
    ///
    /// Returns None if there is no sign change.
    pub fn zero_crossing(&self, node_a: usize, node_b: usize) -> Option<f64> {
        let pa = self.phi[node_a];
        let pb = self.phi[node_b];
        if pa * pb > 0.0 {
            return None;
        }
        let denom = pa - pb;
        if denom.abs() < 1e-60 {
            return Some(0.5);
        }
        Some(pa / denom)
    }
}
/// Types of XFEM enrichment applied to nodes near a crack.
pub enum EnrichmentType {
    /// Heaviside (step) enrichment for nodes whose support is cut by the crack.
    HeavisideJump,
    /// Crack-tip enrichment for nodes near the crack tip.
    CrackTip,
}
/// Registry of phantom nodes for all cracked elements in the mesh.
#[derive(Debug, Clone, Default)]
pub struct PhantomNodeRegistry {
    /// Map from real node index to the corresponding phantom node.
    pub phantoms: Vec<PhantomNode>,
}
impl PhantomNodeRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a phantom node for the given real node.  Returns the phantom index.
    pub fn register(&mut self, real_node_id: usize) -> usize {
        let idx = self.phantoms.len();
        self.phantoms.push(PhantomNode::new(real_node_id));
        idx
    }
    /// Activate the phantom node at `phantom_idx`.
    pub fn activate(&mut self, phantom_idx: usize) {
        if let Some(p) = self.phantoms.get_mut(phantom_idx) {
            p.activate();
        }
    }
    /// Number of active phantom nodes.
    pub fn active_count(&self) -> usize {
        self.phantoms.iter().filter(|p| p.active).count()
    }
    /// Total displacement jump energy (½ ‖Δu‖² summed over active phantoms).
    pub fn jump_energy(&self) -> f64 {
        self.phantoms
            .iter()
            .filter(|p| p.active)
            .map(|p| {
                let j = p.displacement_jump();
                0.5 * (j[0] * j[0] + j[1] * j[1] + j[2] * j[2])
            })
            .sum()
    }
}
/// A 3D crack front represented as an ordered list of points.
///
/// In 3D, the crack front is a curve (not just a point).
#[allow(dead_code)]
pub struct CrackFront3D {
    /// Ordered points along the crack front.
    pub front_points: Vec<[f64; 3]>,
    /// Normal direction to the crack plane at each front point.
    pub normals: Vec<[f64; 3]>,
    /// Tangent direction along the crack front.
    pub tangents: Vec<[f64; 3]>,
}
#[allow(dead_code)]
impl CrackFront3D {
    /// Create a straight crack front along the z-axis.
    pub fn new_straight(start: [f64; 3], end: [f64; 3], n_points: usize) -> Self {
        let n = n_points.max(2);
        let mut front_points = Vec::with_capacity(n);
        let mut normals = Vec::with_capacity(n);
        let mut tangents = Vec::with_capacity(n);
        let dx = end[0] - start[0];
        let dy = end[1] - start[1];
        let dz = end[2] - start[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-30);
        let tang = [dx / len, dy / len, dz / len];
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            front_points.push([start[0] + t * dx, start[1] + t * dy, start[2] + t * dz]);
            normals.push([1.0, 0.0, 0.0]);
            tangents.push(tang);
        }
        Self {
            front_points,
            normals,
            tangents,
        }
    }
    /// Arc length of the crack front.
    pub fn arc_length(&self) -> f64 {
        let mut len = 0.0;
        for i in 1..self.front_points.len() {
            let dx = self.front_points[i][0] - self.front_points[i - 1][0];
            let dy = self.front_points[i][1] - self.front_points[i - 1][1];
            let dz = self.front_points[i][2] - self.front_points[i - 1][2];
            len += (dx * dx + dy * dy + dz * dz).sqrt();
        }
        len
    }
    /// Number of front points.
    pub fn n_points(&self) -> usize {
        self.front_points.len()
    }
    /// Propagate each front point by `da` in the direction of the local normal.
    pub fn propagate_uniform(&mut self, da: f64) {
        for (p, n) in self.front_points.iter_mut().zip(self.normals.iter()) {
            p[0] += da * n[0];
            p[1] += da * n[1];
            p[2] += da * n[2];
        }
    }
}
/// Multiple crack system for interaction analysis.
#[allow(dead_code)]
pub struct MultiCrackSystem {
    /// Vector of cracks.
    pub cracks: Vec<Crack>,
}
#[allow(dead_code)]
impl MultiCrackSystem {
    /// Create a new empty multi-crack system.
    pub fn new() -> Self {
        Self { cracks: Vec::new() }
    }
    /// Add a crack to the system.
    pub fn add_crack(&mut self, crack: Crack) {
        self.cracks.push(crack);
    }
    /// Compute minimum distance between any two crack tips.
    ///
    /// Returns the minimum tip-to-tip distance.
    pub fn min_tip_distance(&self) -> f64 {
        let n = self.cracks.len();
        if n < 2 {
            return f64::INFINITY;
        }
        let mut min_d = f64::MAX;
        for i in 0..n {
            for j in i + 1..n {
                let dx = self.cracks[i].tip[0] - self.cracks[j].tip[0];
                let dy = self.cracks[i].tip[1] - self.cracks[j].tip[1];
                let dz = self.cracks[i].tip[2] - self.cracks[j].tip[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < min_d {
                    min_d = d;
                }
            }
        }
        min_d
    }
    /// Check if cracks are likely to interact (tips within interaction radius).
    pub fn are_interacting(&self, interaction_radius: f64) -> bool {
        self.min_tip_distance() < interaction_radius
    }
    /// Total crack length of all cracks.
    pub fn total_length(&self) -> f64 {
        self.cracks.iter().map(|c| c.length()).sum()
    }
}
/// Utilities for computing Stress Intensity Factors using XFEM.
pub struct XfemSifExtractor {
    /// Young's modulus \[Pa\].
    pub e_modulus: f64,
    /// Poisson's ratio.
    pub poisson: f64,
}
impl XfemSifExtractor {
    /// Create a new SIF extractor.
    pub fn new(e: f64, nu: f64) -> Self {
        Self {
            e_modulus: e,
            poisson: nu,
        }
    }
    /// Extract SIFs using the interaction integral (M-integral) method.
    ///
    /// The interaction integral relates the actual and auxiliary fields:
    /// I = integral_A (sigma_ij^aux * du_i/dx1 + sigma_ij * du_i^aux/dx1
    ///     - W_int * delta_1j) * q_j dA
    ///
    /// where q is a smooth weight function.
    ///
    /// K_I = I^(1) * E' / 2, K_II = I^(2) * E' / 2
    /// where E' = E/(1-nu^2) for plane strain.
    ///
    /// # Arguments
    /// * `interaction_integral_mode1` - Interaction integral for mode I auxiliary field
    /// * `interaction_integral_mode2` - Interaction integral for mode II auxiliary field
    pub fn sif_from_interaction_integral(
        &self,
        interaction_integral_mode1: f64,
        interaction_integral_mode2: f64,
    ) -> (f64, f64) {
        let e_prime = self.e_modulus / (1.0 - self.poisson * self.poisson);
        let ki = interaction_integral_mode1 * e_prime / 2.0;
        let kii = interaction_integral_mode2 * e_prime / 2.0;
        (ki, kii)
    }
    /// Extract SIFs from crack opening displacement (COD).
    ///
    /// Uses the asymptotic relation:
    /// u_y = K_I / mu * sqrt(r / (2 pi)) * (kappa + 1) / 4 * sin(theta/2)
    ///
    /// At theta = pi (upper crack face) and theta = -pi (lower face):
    /// COD = u_y(pi) - u_y(-pi) = K_I * (kappa + 1) / mu * sqrt(r / (2 pi))
    ///
    /// where kappa = 3 - 4*nu for plane strain, mu = E / (2*(1+nu)).
    pub fn ki_from_cod(&self, cod: f64, r: f64) -> f64 {
        let mu = self.e_modulus / (2.0 * (1.0 + self.poisson));
        let kappa = 3.0 - 4.0 * self.poisson;
        let factor = (kappa + 1.0) / mu * (r / (2.0 * PI)).sqrt();
        if factor.abs() < 1e-60 {
            return 0.0;
        }
        cod / factor
    }
    /// Compute the J-integral from interaction integral results.
    ///
    /// J = (K_I^2 + K_II^2) / E' for plane strain.
    pub fn j_from_sif(&self, ki: f64, kii: f64) -> f64 {
        let e_prime = self.e_modulus / (1.0 - self.poisson * self.poisson);
        (ki * ki + kii * kii) / e_prime
    }
}
/// A 2-D / 3-D crack represented as an ordered list of straight segments.
pub struct Crack {
    /// Ordered crack segments from the initial notch to the current tip.
    pub segments: Vec<CrackSegment>,
    /// Current crack-tip position.
    pub tip: [f64; 3],
    /// Direction of most-recent (or initial) growth, unit vector.
    pub growth_direction: [f64; 3],
}
impl Crack {
    /// Create a single-segment crack from `start` to `tip`.
    pub fn new(start: [f64; 3], tip: [f64; 3]) -> Self {
        let dx = tip[0] - start[0];
        let dy = tip[1] - start[1];
        let dz = tip[2] - start[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt().max(f64::EPSILON);
        let growth_direction = [dx / len, dy / len, dz / len];
        let segment = CrackSegment { start, end: tip };
        Crack {
            segments: vec![segment],
            tip,
            growth_direction,
        }
    }
    /// Extend the crack by `da` in the given `direction`.
    pub fn propagate(&mut self, da: f64, direction: [f64; 3]) {
        let dx = direction[0];
        let dy = direction[1];
        let dz = direction[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt().max(f64::EPSILON);
        let unit = [dx / len, dy / len, dz / len];
        let new_tip = [
            self.tip[0] + da * unit[0],
            self.tip[1] + da * unit[1],
            self.tip[2] + da * unit[2],
        ];
        self.segments.push(CrackSegment {
            start: self.tip,
            end: new_tip,
        });
        self.tip = new_tip;
        self.growth_direction = unit;
    }
    /// Total arc-length of the crack.
    pub fn length(&self) -> f64 {
        self.segments.iter().map(|s| s.length()).sum()
    }
    /// Euclidean distance from `point` to the current crack tip.
    pub fn tip_distance(&self, point: [f64; 3]) -> f64 {
        let dx = point[0] - self.tip[0];
        let dy = point[1] - self.tip[1];
        let dz = point[2] - self.tip[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Result of a crack branching check.
#[derive(Debug, Clone, PartialEq)]
pub enum BranchingDecision {
    /// Crack continues along a single propagation path.
    Continue,
    /// Crack bifurcates into two branches at the given angles (radians from tip
    /// direction).
    Bifurcate {
        /// Angle of the first branch (positive = counter-clockwise).
        angle1: f64,
        /// Angle of the second branch.
        angle2: f64,
    },
}
/// Stress intensity factors computed via the modified crack closure integral
/// (MCCI) method.
#[derive(Debug, Clone, Copy)]
pub struct MixedModeSif {
    /// Mode-I stress intensity factor (opening).
    pub k1: f64,
    /// Mode-II stress intensity factor (in-plane shear).
    pub k2: f64,
    /// Mode-III stress intensity factor (out-of-plane shear).
    pub k3: f64,
}
impl MixedModeSif {
    /// Compute effective SIF: K_eff = sqrt(K_I² + K_II² + K_III²/(1-ν))
    pub fn effective(&self, poisson: f64) -> f64 {
        let k3_eff = self.k3 / (1.0 - poisson).max(f64::EPSILON);
        (self.k1 * self.k1 + self.k2 * self.k2 + k3_eff * k3_eff).sqrt()
    }
    /// Compute the energy release rate G = (K_I² + K_II²) / E' + K_III² / (2μ)
    /// where `e_prime` = E / (1 - ν²) and `mu` = E / (2(1+ν)).
    pub fn energy_release_rate(&self, e_prime: f64, shear_modulus: f64) -> f64 {
        (self.k1 * self.k1 + self.k2 * self.k2) / e_prime.max(f64::EPSILON)
            + self.k3 * self.k3 / (2.0 * shear_modulus.max(f64::EPSILON))
    }
}
/// A triangular sub-cell used for integration in cut elements.
///
/// When a crack splits an element, the element is divided into sub-triangles
/// on each side of the crack so that accurate numerical integration can be
/// performed without the integrand discontinuity.
#[derive(Debug, Clone)]
pub struct SubTriangle {
    /// Vertices of the triangle in reference coordinates (xi, eta).
    pub vertices: [[f64; 2]; 3],
    /// Which side of the crack this sub-triangle belongs to (+1 or -1).
    pub side: i8,
}
impl SubTriangle {
    /// Create a new sub-triangle.
    pub fn new(vertices: [[f64; 2]; 3], side: i8) -> Self {
        Self { vertices, side }
    }
    /// Area of the sub-triangle in reference space.
    pub fn area(&self) -> f64 {
        let [v0, v1, v2] = &self.vertices;
        let a = [v1[0] - v0[0], v1[1] - v0[1]];
        let b = [v2[0] - v0[0], v2[1] - v0[1]];
        0.5 * (a[0] * b[1] - a[1] * b[0]).abs()
    }
    /// Centroid of the sub-triangle.
    pub fn centroid(&self) -> [f64; 2] {
        let [v0, v1, v2] = &self.vertices;
        [(v0[0] + v1[0] + v2[0]) / 3.0, (v0[1] + v1[1] + v2[1]) / 3.0]
    }
    /// 3-point Gauss quadrature points and weights in barycentric coordinates,
    /// mapped back to (xi, eta) reference space.
    ///
    /// Returns a list of `(xi, eta, weight)` tuples.
    pub fn gauss_points(&self) -> [(f64, f64, f64); 3] {
        let area = self.area();
        let w = area / 3.0;
        let [v0, v1, v2] = &self.vertices;
        let pts = [
            [(v0[0] + v1[0]) / 2.0, (v0[1] + v1[1]) / 2.0],
            [(v1[0] + v2[0]) / 2.0, (v1[1] + v2[1]) / 2.0],
            [(v2[0] + v0[0]) / 2.0, (v2[1] + v0[1]) / 2.0],
        ];
        [
            (pts[0][0], pts[0][1], w),
            (pts[1][0], pts[1][1], w),
            (pts[2][0], pts[2][1], w),
        ]
    }
}
/// A node along a 3-D crack front that carries branch-function enrichment DOFs.
#[derive(Debug, Clone)]
pub struct CrackFrontNode {
    /// 3-D position of this front node.
    pub position: [f64; 3],
    /// Unit tangent to the crack front at this node.
    pub tangent: [f64; 3],
    /// Branch function enrichment DOFs: four per node for the standard
    /// Moës–Dolbow–Belytschko tip enrichment {√r sinθ/2, √r cosθ/2,
    ///  √r sinθ/2 sinθ, √r cosθ/2 sinθ}.
    pub branch_dofs: [f64; 4],
    /// Local crack-front coordinate (arc length from the start of the front).
    pub arc_length: f64,
}
impl CrackFrontNode {
    /// Create a new front node with zero enrichment DOFs.
    pub fn new(position: [f64; 3], tangent: [f64; 3], arc_length: f64) -> Self {
        let t_norm = vec3_norm(tangent);
        let t_unit = if t_norm > f64::EPSILON {
            [
                tangent[0] / t_norm,
                tangent[1] / t_norm,
                tangent[2] / t_norm,
            ]
        } else {
            [1.0, 0.0, 0.0]
        };
        Self {
            position,
            tangent: t_unit,
            branch_dofs: [0.0; 4],
            arc_length,
        }
    }
    /// Evaluate the four branch (tip enrichment) functions at local polar
    /// coordinates `(r, theta)` around this front node.
    pub fn branch_functions(r: f64, theta: f64) -> [f64; 4] {
        let sr = r.sqrt();
        let ht = theta / 2.0;
        [
            sr * ht.sin(),
            sr * ht.cos(),
            sr * ht.sin() * theta.sin(),
            sr * ht.cos() * theta.sin(),
        ]
    }
    /// Enriched displacement contribution from this front node at `(r, theta)`.
    pub fn enriched_displacement(&self, r: f64, theta: f64) -> [f64; 4] {
        let bf = Self::branch_functions(r, theta);
        [
            self.branch_dofs[0] * bf[0],
            self.branch_dofs[1] * bf[1],
            self.branch_dofs[2] * bf[2],
            self.branch_dofs[3] * bf[3],
        ]
    }
}
/// A discrete 3-D crack front represented as an ordered list of enrichment nodes.
///
/// This is distinct from [`CrackFront3D`] (which stores raw geometry points)
/// and focuses on XFEM branch-function enrichment.
#[derive(Debug, Clone, Default)]
pub struct CrackFrontEnrichment {
    /// Ordered list of nodes along the crack front.
    pub nodes: Vec<CrackFrontNode>,
}
impl CrackFrontEnrichment {
    /// Create an empty crack front.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a new front node, automatically computing the arc-length increment.
    pub fn push_node(&mut self, position: [f64; 3], tangent: [f64; 3]) {
        let arc = if let Some(last) = self.nodes.last() {
            let d = [
                position[0] - last.position[0],
                position[1] - last.position[1],
                position[2] - last.position[2],
            ];
            last.arc_length + vec3_norm(d)
        } else {
            0.0
        };
        self.nodes.push(CrackFrontNode::new(position, tangent, arc));
    }
    /// Total arc length of the front.
    pub fn total_arc_length(&self) -> f64 {
        self.nodes.last().map(|n| n.arc_length).unwrap_or(0.0)
    }
    /// Number of nodes along the front.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }
    /// Returns true if the front has no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
    /// Interpolate the crack-front position at arc-length parameter `s` using
    /// linear interpolation between adjacent nodes.
    pub fn position_at(&self, s: f64) -> Option<[f64; 3]> {
        if self.nodes.is_empty() {
            return None;
        }
        let s = s.clamp(0.0, self.total_arc_length());
        for i in 1..self.nodes.len() {
            let n0 = &self.nodes[i - 1];
            let n1 = &self.nodes[i];
            if s <= n1.arc_length {
                let t = (s - n0.arc_length) / (n1.arc_length - n0.arc_length).max(f64::EPSILON);
                return Some([
                    n0.position[0] + t * (n1.position[0] - n0.position[0]),
                    n0.position[1] + t * (n1.position[1] - n0.position[1]),
                    n0.position[2] + t * (n1.position[2] - n0.position[2]),
                ]);
            }
        }
        self.nodes.last().map(|n| n.position)
    }
}
/// A single straight crack segment defined by two endpoints in 3-D space.
pub struct CrackSegment {
    /// Start point of the segment.
    pub start: [f64; 3],
    /// End point of the segment.
    pub end: [f64; 3],
}
impl CrackSegment {
    /// Euclidean length of this segment.
    pub fn length(&self) -> f64 {
        let dx = self.end[0] - self.start[0];
        let dy = self.end[1] - self.start[1];
        let dz = self.end[2] - self.start[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// A phantom node used to model displacement discontinuities across crack faces.
///
/// In the phantom node method a cracked element is replaced by two overlapping
/// sub-elements, each carrying one half of the original element plus a phantom
/// (fictitious) counterpart.  The phantom nodes hold the displacement jump DOFs.
#[derive(Debug, Clone)]
pub struct PhantomNode {
    /// Index of the real node that this phantom duplicates.
    pub real_node_id: usize,
    /// Phantom DOF values (displacement enrichment).
    pub dofs: [f64; 3],
    /// Whether this phantom node is active (crack has passed through it).
    pub active: bool,
}
impl PhantomNode {
    /// Create a new inactive phantom node mirroring `real_node_id`.
    pub fn new(real_node_id: usize) -> Self {
        Self {
            real_node_id,
            dofs: [0.0; 3],
            active: false,
        }
    }
    /// Activate this phantom node (called when crack front enters the element).
    pub fn activate(&mut self) {
        self.active = true;
    }
    /// Compute the displacement jump `u+ - u-` at this node.
    pub fn displacement_jump(&self) -> [f64; 3] {
        [self.dofs[0] * 2.0, self.dofs[1] * 2.0, self.dofs[2] * 2.0]
    }
}
