// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Coupling method for meshfree-FEM interface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CouplingMethod {
    /// Direct coupling: meshfree nodes on FE boundary are tied.
    Direct,
    /// Bridging domain method with overlap region.
    BridgingDomain,
    /// Arlequin method with energy partition.
    Arlequin,
}
/// RKPM shape function result.
#[derive(Debug, Clone)]
pub struct RkpmShape {
    /// Shape function values.
    pub psi: Vec<f64>,
    /// d(psi)/dx.
    pub dpsi_dx: Vec<f64>,
    /// d(psi)/dy.
    pub dpsi_dy: Vec<f64>,
    /// Neighbour indices.
    pub neighbours: Vec<usize>,
}
/// Result of a convergence rate estimation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConvergenceRate {
    /// Estimated order of convergence for L2 norm.
    pub l2_rate: f64,
    /// Estimated order of convergence for L-infinity norm.
    pub linf_rate: f64,
    /// Estimated order of convergence for H1 semi-norm.
    pub h1_rate: f64,
    /// Estimated order of convergence for energy norm.
    pub energy_rate: f64,
}
/// Recovered stress at a point (plane stress: sigma_xx, sigma_yy, tau_xy).
#[derive(Debug, Clone, Copy, Default)]
pub struct StressState {
    /// Normal stress in x.
    pub sigma_xx: f64,
    /// Normal stress in y.
    pub sigma_yy: f64,
    /// Shear stress.
    pub tau_xy: f64,
}
impl StressState {
    /// Von Mises equivalent stress for plane stress.
    #[allow(dead_code)]
    pub fn von_mises(&self) -> f64 {
        (self.sigma_xx * self.sigma_xx + self.sigma_yy * self.sigma_yy
            - self.sigma_xx * self.sigma_yy
            + 3.0 * self.tau_xy * self.tau_xy)
            .sqrt()
    }
    /// Principal stresses (max, min).
    #[allow(dead_code)]
    pub fn principal(&self) -> (f64, f64) {
        let avg = (self.sigma_xx + self.sigma_yy) / 2.0;
        let diff = (self.sigma_xx - self.sigma_yy) / 2.0;
        let r = (diff * diff + self.tau_xy * self.tau_xy).sqrt();
        (avg + r, avg - r)
    }
}
/// RBF interpolant for scattered 2-D data.
#[derive(Debug, Clone)]
pub struct RbfInterpolant {
    /// Centre positions.
    pub centres: Vec<Point2>,
    /// Coefficients (one per centre + polynomial augmentation).
    pub coeffs: Vec<f64>,
    /// Shape parameter.
    pub c: f64,
    /// RBF type.
    pub rtype: RbfType,
}
impl RbfInterpolant {
    /// Evaluate the interpolant at a given point.
    #[allow(dead_code)]
    pub fn eval(&self, pt: Point2) -> f64 {
        let mut val = 0.0;
        for (i, centre) in self.centres.iter().enumerate() {
            let r = dist2(pt, *centre);
            val += self.coeffs[i] * rbf_eval(r, self.c, self.rtype);
        }
        val
    }
}
/// Summary of a convergence study for display / reporting.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConvergenceStudySummary {
    /// All data points in the study.
    pub points: Vec<ConvergencePoint>,
    /// Computed rates between successive refinements.
    pub rates: Vec<ConvergenceRate>,
    /// Average L2 convergence rate.
    pub avg_l2_rate: f64,
    /// Average H1 convergence rate.
    pub avg_h1_rate: f64,
}
/// Coupling interface between meshfree and FEM domains.
#[derive(Debug, Clone)]
pub struct CouplingInterface {
    /// Interface nodes.
    pub nodes: Vec<CouplingNode>,
    /// Coupling method.
    pub method: CouplingMethod,
    /// Penalty parameter for direct coupling.
    pub penalty: f64,
}
impl CouplingInterface {
    /// Create a new coupling interface.
    #[allow(dead_code)]
    pub fn new(method: CouplingMethod, penalty: f64) -> Self {
        Self {
            nodes: Vec::new(),
            method,
            penalty,
        }
    }
    /// Add a coupling node.
    #[allow(dead_code)]
    pub fn add_node(&mut self, node: CouplingNode) {
        self.nodes.push(node);
    }
    /// Compute the coupling constraint matrix for direct coupling.
    ///
    /// Returns a penalty contribution to be added to the combined stiffness.
    #[allow(dead_code)]
    pub fn compute_penalty_coupling(&self, ndof_meshfree: usize, ndof_fem: usize) -> Vec<f64> {
        let ndof_total = ndof_meshfree + ndof_fem;
        let mut coupling_k = vec![0.0; ndof_total * ndof_total];
        for cn in &self.nodes {
            if let (Some(mi), Some(fi)) = (cn.meshfree_idx, cn.fem_idx) {
                for dof in 0..2 {
                    let im = mi * 2 + dof;
                    let jf = ndof_meshfree + fi * 2 + dof;
                    if im < ndof_meshfree && jf < ndof_total {
                        coupling_k[im * ndof_total + im] += self.penalty;
                        coupling_k[jf * ndof_total + jf] += self.penalty;
                        coupling_k[im * ndof_total + jf] -= self.penalty;
                        coupling_k[jf * ndof_total + im] -= self.penalty;
                    }
                }
            }
        }
        coupling_k
    }
    /// Compute blending weights for the bridging domain method.
    ///
    /// In the overlap region, the total energy is partitioned:
    ///   E_total = w * E_meshfree + (1-w) * E_fem
    #[allow(dead_code)]
    pub fn compute_bridging_weights(&mut self, overlap_width: f64) {
        for cn in &mut self.nodes {
            cn.blend_weight = cn.blend_weight.clamp(0.0, 1.0);
            let _ = overlap_width;
        }
    }
    /// Number of coupling nodes.
    #[allow(dead_code)]
    pub fn num_coupling_nodes(&self) -> usize {
        self.nodes.len()
    }
}
/// Support domain manager — brute-force neighbour search.
///
/// For large problems a cell-list or kd-tree should replace this; the brute
/// force version is kept for clarity and correctness.
#[derive(Debug, Clone)]
pub struct SupportDomain {
    /// All nodes in the cloud.
    pub nodes: Vec<MeshfreeNode>,
    /// Global dilation parameter (multiplied with each node's radius).
    pub dilation: f64,
}
impl SupportDomain {
    /// Create a new support domain from a set of nodes and a dilation factor.
    #[allow(dead_code)]
    pub fn new(nodes: Vec<MeshfreeNode>, dilation: f64) -> Self {
        Self { nodes, dilation }
    }
    /// Find indices of all nodes whose support covers `pt`.
    #[allow(dead_code)]
    pub fn neighbours(&self, pt: Point2) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                let d = dist2(pt, n.pos);
                if d < n.support_radius * self.dilation {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Find indices of all nodes whose support covers `pt`, returning
    /// `(index, distance)` pairs sorted by distance.
    #[allow(dead_code)]
    pub fn neighbours_sorted(&self, pt: Point2) -> Vec<(usize, f64)> {
        let mut pairs: Vec<(usize, f64)> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                let d = dist2(pt, n.pos);
                if d < n.support_radius * self.dilation {
                    Some((i, d))
                } else {
                    None
                }
            })
            .collect();
        pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        pairs
    }
    /// Return the number of nodes in the cloud.
    #[allow(dead_code)]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}
/// Result from an EFG plane-stress analysis.
#[derive(Debug, Clone)]
pub struct EfgResult {
    /// Nodal displacements (2 DOF per node: ux, uy).
    pub displacements: Vec<f64>,
    /// Number of nodes.
    pub num_nodes: usize,
}
/// A node in a meshfree particle cloud.
#[derive(Debug, Clone)]
pub struct MeshfreeNode {
    /// Position in 2-D.
    pub pos: Point2,
    /// Support radius (influence domain size).
    pub support_radius: f64,
    /// Nodal value (scalar field).
    pub value: f64,
}
/// Weight (kernel) function type for meshfree approximation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelType {
    /// Cubic B-spline kernel.
    CubicSpline,
    /// Quartic spline kernel.
    QuarticSpline,
    /// Gaussian kernel with shape parameter.
    Gaussian,
    /// Wendland C2 compactly-supported kernel.
    WendlandC2,
    /// Wendland C4 compactly-supported kernel.
    WendlandC4,
}
/// An essential boundary condition specification.
#[derive(Debug, Clone)]
pub struct EssentialBc {
    /// Node index.
    pub node: usize,
    /// DOF index (0 = x, 1 = y).
    pub dof: usize,
    /// Prescribed value.
    pub value: f64,
}
/// Method for enforcing essential BCs in meshfree methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcEnforcementMethod {
    /// Penalty method.
    Penalty,
    /// Lagrange multiplier method.
    LagrangeMultiplier,
    /// Nitsche's method (penalty-free weak enforcement).
    Nitsche,
}
/// RBF type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RbfType {
    /// Multiquadric: sqrt(r^2 + c^2).
    Multiquadric,
    /// Inverse multiquadric: 1 / sqrt(r^2 + c^2).
    InverseMultiquadric,
    /// Gaussian: exp(-c^2 r^2).
    Gaussian,
    /// Thin-plate spline: r^2 ln(r).
    ThinPlateSpline,
    /// Polyharmonic spline r^3.
    Polyharmonic3,
}
/// A finite element represented by its node positions and connectivity.
#[derive(Debug, Clone)]
pub struct FemElement {
    /// Node positions (2-D).
    pub nodes: Vec<Point2>,
    /// Element connectivity (indices into the global node list).
    pub connectivity: Vec<usize>,
    /// Element stiffness (dense, row-major).
    pub stiffness: Vec<f64>,
}
/// Error indicator for a single node.
#[derive(Debug, Clone, Copy)]
pub struct ErrorIndicator {
    /// Node index.
    pub node: usize,
    /// Estimated error (e.g. von Mises stress gradient magnitude).
    pub error: f64,
}
/// Interface node for meshfree-FEM coupling.
#[derive(Debug, Clone)]
pub struct CouplingNode {
    /// Position.
    pub pos: Point2,
    /// Meshfree node index (if exists).
    pub meshfree_idx: Option<usize>,
    /// FEM node index (if exists).
    pub fem_idx: Option<usize>,
    /// Weight for blending (0 = pure FEM, 1 = pure meshfree).
    pub blend_weight: f64,
}
/// A triangular background integration cell.
#[derive(Debug, Clone, Copy)]
pub struct TriCell {
    /// Vertex indices (into integration point list or node list).
    pub vertices: [Point2; 3],
}
impl TriCell {
    /// Area of the triangle.
    #[allow(dead_code)]
    pub fn area(&self) -> f64 {
        let v = self.vertices;
        0.5 * ((v[1][0] - v[0][0]) * (v[2][1] - v[0][1])
            - (v[2][0] - v[0][0]) * (v[1][1] - v[0][1]))
            .abs()
    }
    /// Centroid of the triangle.
    #[allow(dead_code)]
    pub fn centroid(&self) -> Point2 {
        let v = self.vertices;
        [
            (v[0][0] + v[1][0] + v[2][0]) / 3.0,
            (v[0][1] + v[1][1] + v[2][1]) / 3.0,
        ]
    }
}
/// Material parameters for plane-stress EFG analysis.
#[derive(Debug, Clone, Copy)]
pub struct EfgMaterial {
    /// Young's modulus.
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Thickness (for plane stress).
    pub thickness: f64,
}
/// A single data point in a convergence study.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConvergencePoint {
    /// Characteristic node spacing h.
    pub h: f64,
    /// Number of nodes in the discretisation.
    pub num_nodes: usize,
    /// L2 error norm.
    pub l2_error: f64,
    /// L-infinity (max) error norm.
    pub linf_error: f64,
    /// H1 semi-norm error (gradient error).
    pub h1_error: f64,
    /// Energy norm error.
    pub energy_error: f64,
}
/// Result of MLS shape function evaluation at a single point.
#[derive(Debug, Clone)]
pub struct MlsShape {
    /// Shape function values — one per neighbour.
    pub phi: Vec<f64>,
    /// d(phi)/dx for each neighbour.
    pub dphi_dx: Vec<f64>,
    /// d(phi)/dy for each neighbour.
    pub dphi_dy: Vec<f64>,
    /// Neighbour indices.
    pub neighbours: Vec<usize>,
}
