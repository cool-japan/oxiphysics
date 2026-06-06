//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute the penalty contact force on a node that has penetrated a half-space.
pub fn nodal_contact_force(
    node_pos: [f64; 3],
    obstacle_normal: [f64; 3],
    obstacle_d: f64,
    penalty: f64,
) -> [f64; 3] {
    let gap = gap_function_node_to_plane(node_pos, obstacle_normal, obstacle_d);
    let f_mag = penalty * (-gap).max(0.0);
    [
        f_mag * obstacle_normal[0],
        f_mag * obstacle_normal[1],
        f_mag * obstacle_normal[2],
    ]
}
/// Compute the signed gap from a node to a plane.
pub fn gap_function_node_to_plane(node: [f64; 3], plane_normal: [f64; 3], plane_d: f64) -> f64 {
    plane_normal[0] * node[0] + plane_normal[1] * node[1] + plane_normal[2] * node[2] + plane_d
}
/// Compute the signed gap and contact normal from a node to the nearest point on a line segment.
pub fn gap_function_node_to_segment(node: [f64; 3], p0: [f64; 3], p1: [f64; 3]) -> (f64, [f64; 3]) {
    let ab = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let ap = [node[0] - p0[0], node[1] - p0[1], node[2] - p0[2]];
    let ab_len_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let t = if ab_len_sq > 0.0 {
        let dot = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2];
        (dot / ab_len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let nearest = [p0[0] + t * ab[0], p0[1] + t * ab[1], p0[2] + t * ab[2]];
    let diff = [
        node[0] - nearest[0],
        node[1] - nearest[1],
        node[2] - nearest[2],
    ];
    let gap = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
    let normal = if gap > 0.0 {
        [diff[0] / gap, diff[1] / gap, diff[2] / gap]
    } else {
        [0.0, 0.0, 0.0]
    };
    (gap, normal)
}
/// Find candidate self-contact node pairs within a single mesh.
///
/// Returns pairs `(i, j)` where `i < j` and the nodes are within `cutoff`.
/// Nodes are not paired with themselves.
pub fn self_contact_candidate_pairs(nodes: &[[f64; 3]], cutoff: f64) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let n = nodes.len();
    for i in 0..n {
        for j in i + 1..n {
            let dx = nodes[i][0] - nodes[j][0];
            let dy = nodes[i][1] - nodes[j][1];
            let dz = nodes[i][2] - nodes[j][2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist <= cutoff {
                pairs.push((i, j));
            }
        }
    }
    pairs
}
/// Compute the gap and contact normal for a node versus a triangular surface element.
///
/// Projects the node onto the triangle plane and computes the signed distance
/// (positive = above the triangle on the normal side).
/// The normal is the outward unit normal of the triangle.
///
/// Returns `(gap, normal)` where `gap > 0` means the node is on the outside.
pub fn node_to_triangle_gap(
    node: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> (f64, [f64; 3]) {
    let e0 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e1 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    let n = [
        e0[1] * e1[2] - e0[2] * e1[1],
        e0[2] * e1[0] - e0[0] * e1[2],
        e0[0] * e1[1] - e0[1] * e1[0],
    ];
    let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if n_len < 1e-30 {
        return (0.0, [0.0, 0.0, 0.0]);
    }
    let normal = [n[0] / n_len, n[1] / n_len, n[2] / n_len];
    let pv = [node[0] - v0[0], node[1] - v0[1], node[2] - v0[2]];
    let gap = pv[0] * normal[0] + pv[1] * normal[1] + pv[2] * normal[2];
    (gap, normal)
}
/// Compute the area of a triangle from three vertices.
pub fn triangle_area(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let e0 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e1 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    let cross = [
        e0[1] * e1[2] - e0[2] * e1[1],
        e0[2] * e1[0] - e0[0] * e1[2],
        e0[0] * e1[1] - e0[1] * e1[0],
    ];
    0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
}
/// Compute the 3×3 penalty contact stiffness matrix for a single node-to-plane
/// contact pair: `K = k_n * n ⊗ n`.
pub fn penalty_contact_stiffness(k_n: f64, normal: [f64; 3]) -> [[f64; 3]; 3] {
    let mut ke = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            ke[i][j] = k_n * normal[i] * normal[j];
        }
    }
    ke
}
/// Estimate the contact patch area from nodal gaps and element area.
///
/// Counts nodes in contact (gap < 0) and multiplies by the average element
/// area to estimate the contact patch area.
pub fn contact_patch_area(gaps: &[f64], element_area: f64) -> f64 {
    let n_contact = gaps.iter().filter(|&&g| g < 0.0).count();
    n_contact as f64 * element_area
}
/// Hertz contact radius for two elastic spheres.
///
/// `a = (3 * F * R* / (4 * E*))^(1/3)`
///
/// where `R* = (1/R1 + 1/R2)^-1` (here `R1 = R2 = r`) and
/// `E* = ((1-nu1^2)/E1 + (1-nu2^2)/E2)^-1`.
pub fn hertz_contact_radius(force: f64, radius: f64, e1: f64, e2: f64, nu1: f64, nu2: f64) -> f64 {
    let r_star = radius / 2.0;
    let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
    (3.0 * force * r_star / (4.0 * e_star)).powf(1.0 / 3.0)
}
/// Hertz maximum contact pressure: `p0 = 3F / (2 * pi * a^2)`.
pub fn hertz_max_pressure(force: f64, radius: f64, e1: f64, e2: f64, nu1: f64, nu2: f64) -> f64 {
    let a = hertz_contact_radius(force, radius, e1, e2, nu1, nu2);
    if a < 1e-30 {
        return 0.0;
    }
    3.0 * force / (2.0 * std::f64::consts::PI * a * a)
}
/// Hertz contact indentation depth (approach): `delta = a^2 / R*`.
pub fn hertz_indentation(force: f64, radius: f64, e1: f64, e2: f64, nu1: f64, nu2: f64) -> f64 {
    let r_star = radius / 2.0;
    let a = hertz_contact_radius(force, radius, e1, e2, nu1, nu2);
    a * a / r_star
}
/// Update gap function after large deformation using current node position.
///
/// The deformed gap accounts for rigid-body rotation of the contact surface.
pub fn deformed_gap_plane(
    node_pos: [f64; 3],
    deformed_normal: [f64; 3],
    plane_point: [f64; 3],
) -> f64 {
    let dp = [
        node_pos[0] - plane_point[0],
        node_pos[1] - plane_point[1],
        node_pos[2] - plane_point[2],
    ];
    dp[0] * deformed_normal[0] + dp[1] * deformed_normal[1] + dp[2] * deformed_normal[2]
}
/// Compute the incremental contact force for a large-deformation step.
///
/// Accounts for the change in normal direction between timesteps.
pub fn large_deformation_contact_force(
    gap: f64,
    normal_old: [f64; 3],
    normal_new: [f64; 3],
    penalty: f64,
    lambda: f64,
) -> [f64; 3] {
    let pen = (-gap).max(0.0);
    let f_mag = lambda + penalty * pen;
    let n = if gap < 0.0 { normal_new } else { normal_old };
    [f_mag * n[0], f_mag * n[1], f_mag * n[2]]
}
/// Project a point onto a line segment \[a, b\], returning the projected point.
pub fn project_point_to_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let ab_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    if ab_sq < 1e-30 {
        return a;
    }
    let t = ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / ab_sq).clamp(0.0, 1.0);
    [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]]
}
/// 1D Gaussian quadrature points and weights on \[-1, 1\].
pub(super) fn gauss_quad_1d(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![0.0], vec![2.0]),
        2 => (
            vec![-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.7745966692, 0.0, 0.7745966692],
            vec![0.5555555556, 0.8888888889, 0.5555555556],
        ),
        4 => (
            vec![-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116],
            vec![0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451],
        ),
        _ => {
            let pts: Vec<f64> = (0..n)
                .map(|i| -1.0 + (2.0 * i as f64 + 1.0) / n as f64)
                .collect();
            let w = 2.0 / n as f64;
            let wts = vec![w; n];
            (pts, wts)
        }
    }
}
/// Gap function for node against a quadrilateral face (4-node bilinear).
///
/// Returns signed distance from node to bilinear surface.
pub fn gap_function_node_to_quad(
    node: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
) -> (f64, [f64; 3]) {
    let (g1, n1) = node_to_triangle_gap_fem(node, v0, v1, v2);
    let (g2, n2) = node_to_triangle_gap_fem(node, v0, v2, v3);
    let g = 0.5 * (g1 + g2);
    let n = [
        0.5 * (n1[0] + n2[0]),
        0.5 * (n1[1] + n2[1]),
        0.5 * (n1[2] + n2[2]),
    ];
    (g, n)
}
/// Compute gap from node to triangle (signed distance to triangle plane).
pub(super) fn node_to_triangle_gap_fem(
    node: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> (f64, [f64; 3]) {
    let e0 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e1 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    let n = [
        e0[1] * e1[2] - e0[2] * e1[1],
        e0[2] * e1[0] - e0[0] * e1[2],
        e0[0] * e1[1] - e0[1] * e1[0],
    ];
    let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if n_len < 1e-30 {
        return (0.0, [0.0, 0.0, 0.0]);
    }
    let normal = [n[0] / n_len, n[1] / n_len, n[2] / n_len];
    let pv = [node[0] - v0[0], node[1] - v0[1], node[2] - v0[2]];
    let gap = pv[0] * normal[0] + pv[1] * normal[1] + pv[2] * normal[2];
    (gap, normal)
}
/// Compute contact traction (force per unit area) from penalty parameters.
pub fn contact_traction(
    gap: f64,
    tangential_disp: f64,
    penalty_normal: f64,
    penalty_tangential: f64,
    friction_coeff: f64,
) -> (f64, f64) {
    let t_n = penalty_normal * (-gap).max(0.0);
    let t_t_trial = penalty_tangential * tangential_disp;
    let limit = friction_coeff * t_n;
    let t_t = if t_t_trial.abs() <= limit {
        t_t_trial
    } else {
        limit * t_t_trial.signum()
    };
    (t_n, t_t)
}
/// Reconstruct contact pressure distribution from nodal contact forces.
///
/// Uses a simple inverse-distance weighting to smooth nodal forces
/// onto a spatial pressure field.
pub fn contact_pressure_field(
    query_points: &[[f64; 3]],
    contact_nodes: &[[f64; 3]],
    nodal_forces: &[f64],
    element_area: f64,
) -> Vec<f64> {
    if contact_nodes.is_empty() || element_area < 1e-30 {
        return vec![0.0; query_points.len()];
    }
    let nodal_pressures: Vec<f64> = nodal_forces.iter().map(|&f| f / element_area).collect();
    query_points
        .iter()
        .map(|&q| {
            let mut sum_w = 0.0;
            let mut sum_wp = 0.0;
            for (node, &p) in contact_nodes.iter().zip(nodal_pressures.iter()) {
                let dx = q[0] - node[0];
                let dy = q[1] - node[1];
                let dz = q[2] - node[2];
                let dist2 = dx * dx + dy * dy + dz * dz;
                let w = if dist2 < 1e-30 { 1e30 } else { 1.0 / dist2 };
                sum_w += w;
                sum_wp += w * p;
            }
            if sum_w < 1e-30 { 0.0 } else { sum_wp / sum_w }
        })
        .collect()
}
/// Assemble the contact contribution into a global stiffness matrix.
///
/// For each contact pair, adds the penalty stiffness contribution
/// to the relevant diagonal entries.
///
/// `n_dofs` × `n_dofs` sparse matrix (stored as flat Vec).
pub fn assemble_contact_stiffness(
    n_dofs: usize,
    contact_pairs: &[(usize, usize)],
    gaps: &[f64],
    penalty: f64,
) -> Vec<f64> {
    let mut k = vec![0.0f64; n_dofs * n_dofs];
    for (&(i, j), &gap) in contact_pairs.iter().zip(gaps.iter()) {
        if gap >= 0.0 {
            continue;
        }
        if i < n_dofs {
            k[i * n_dofs + i] += penalty;
        }
        if j < n_dofs {
            k[j * n_dofs + j] += penalty;
        }
        if i < n_dofs && j < n_dofs {
            k[i * n_dofs + j] -= penalty;
            k[j * n_dofs + i] -= penalty;
        }
    }
    k
}
/// Assemble the contact force vector from penalty contact.
pub fn assemble_contact_force(
    n_dofs: usize,
    contact_nodes: &[usize],
    gaps: &[f64],
    penalty: f64,
) -> Vec<f64> {
    let mut f = vec![0.0f64; n_dofs];
    for (&node, &gap) in contact_nodes.iter().zip(gaps.iter()) {
        if gap < 0.0 && node < n_dofs {
            f[node] += penalty * (-gap);
        }
    }
    f
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::contact_fem::*;
    #[test]
    fn test_penalty_contact_normal_force() {
        let pc = PenaltyContact {
            stiffness: 1000.0,
            friction_coeff: 0.3,
        };
        assert_eq!(pc.normal_force(0.05), 0.0);
        let gap = -0.01_f64;
        let f_n = pc.normal_force(gap);
        assert!(f_n > 0.0);
        let expected = 1000.0 * 0.01;
        assert!((f_n - expected).abs() < 1e-12);
    }
    #[test]
    fn test_augmented_lagrangian_multiplier_update() {
        let mut alc = AugmentedLagrangianContact::new(500.0);
        assert_eq!(alc.lambda_n, 0.0);
        alc.update_multipliers(-0.02, 0.0, 0.3);
        assert!(alc.lambda_n < 0.0);
        let expected = 500.0 * (-0.02_f64);
        assert!((alc.lambda_n - expected).abs() < 1e-12);
    }
    #[test]
    fn test_nodal_contact_force_below_plane() {
        let node_pos = [0.0, 0.0, -0.01];
        let normal = [0.0, 0.0, 1.0];
        let d = 0.0;
        let penalty = 1000.0;
        let force = nodal_contact_force(node_pos, normal, d, penalty);
        assert!(force[2] > 0.0);
        assert_eq!(force[0], 0.0);
        assert_eq!(force[1], 0.0);
        let expected_fz = penalty * 0.01;
        assert!((force[2] - expected_fz).abs() < 1e-12);
    }
    #[test]
    fn test_contact_search_brute_force() {
        let nodes_a: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let nodes_b: &[[f64; 3]] = &[[0.1, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let cutoff = 0.5;
        let search = ContactSearch::brute_force_search(nodes_a, nodes_b, cutoff);
        assert_eq!(search.pair_count(), 1);
        assert_eq!(search.candidate_pairs[0], (0, 0));
    }
    #[test]
    fn test_gap_function_node_on_plane() {
        let node = [1.0, 2.0, 0.0];
        let normal = [0.0, 0.0, 1.0];
        let gap = gap_function_node_to_plane(node, normal, 0.0);
        assert!((gap).abs() < 1e-15);
    }
    #[test]
    fn test_gap_function_node_to_segment() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [2.0, 0.0, 0.0];
        let node = [1.0, 0.0, 1.0];
        let (gap, normal) = gap_function_node_to_segment(node, p0, p1);
        assert!((gap - 1.0).abs() < 1e-12);
        assert!((normal[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_contact_residual_separated() {
        let pc = PenaltyContact {
            stiffness: 1000.0,
            friction_coeff: 0.3,
        };
        let (f_n, f_t) = pc.contact_residual(0.1, 0.5);
        assert_eq!(f_n, 0.0);
        assert_eq!(f_t, 0.0);
    }
    #[test]
    fn test_augmented_lagrangian_force() {
        let mut alc = AugmentedLagrangianContact::new(200.0);
        let gap = -0.05;
        alc.lambda_n = -3.0;
        let f = alc.force(gap);
        let expected = -3.0 + 200.0 * (-0.05_f64);
        assert!((f - expected).abs() < 1e-12);
    }
    #[test]
    fn test_tied_contact_force() {
        let pairs = vec![(0, 1)];
        let tc = TiedContact::new(pairs, 1000.0);
        let slave_disp = [0.01, 0.0, 0.0];
        let master_disp = [0.0, 0.0, 0.0];
        let (f_s, f_m) = tc.tied_force(slave_disp, master_disp);
        assert!((f_s[0] - (-10.0)).abs() < 1e-12);
        assert!((f_m[0] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_tied_contact_zero_gap() {
        let pairs = vec![(0, 1)];
        let tc = TiedContact::new(pairs, 1000.0);
        let disp = [0.05, 0.02, 0.01];
        let (f_s, f_m) = tc.tied_force(disp, disp);
        for (i, (&fs_v, &fm_v)) in f_s.iter().zip(f_m.iter()).enumerate() {
            assert_eq!(fs_v, 0.0, "zero gap should give zero force at [{i}]");
            assert_eq!(fm_v, 0.0, "zero gap master force should be zero at [{i}]");
        }
    }
    #[test]
    fn test_tied_stiffness_symmetry() {
        let tc = TiedContact::new(vec![], 1000.0);
        let ke = tc.tied_stiffness();
        for (i, row) in ke.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert!(
                    (v - ke[j][i]).abs() < 1e-12,
                    "stiffness should be symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_sliding_contact_stick() {
        let sc = SlidingContact::new(1000.0, 500.0, 0.3, 1);
        let gap = -0.01;
        let normal = [0.0, 0.0, 1.0];
        let tang_disp = [1e-8, 0.0, 0.0];
        let (f_n, f_t, is_sliding) = sc.contact_forces(gap, normal, tang_disp);
        assert!(
            !is_sliding,
            "should be in stick for small tangential displacement"
        );
        assert!(f_n[2] > 0.0, "normal force should be in normal direction");
        assert!(f_t[0].abs() > 0.0, "should have tangential force in stick");
    }
    #[test]
    fn test_sliding_contact_slip() {
        let sc = SlidingContact::new(1000.0, 500.0, 0.3, 1);
        let gap = -0.01;
        let normal = [0.0, 0.0, 1.0];
        let tang_disp = [1.0, 0.0, 0.0];
        let (_f_n, f_t, is_sliding) = sc.contact_forces(gap, normal, tang_disp);
        assert!(
            is_sliding,
            "should be sliding for large tangential displacement"
        );
        let f_n_mag = 1000.0 * 0.01;
        let f_t_mag = (f_t[0] * f_t[0] + f_t[1] * f_t[1] + f_t[2] * f_t[2]).sqrt();
        let limit = 0.3 * f_n_mag;
        assert!(
            (f_t_mag - limit).abs() / limit < 1e-10,
            "f_t_mag = {f_t_mag}, limit = {limit}"
        );
    }
    #[test]
    fn test_sliding_contact_separated() {
        let sc = SlidingContact::new(1000.0, 500.0, 0.3, 1);
        let gap = 0.01;
        let normal = [0.0, 0.0, 1.0];
        let tang_disp = [0.1, 0.0, 0.0];
        let (f_n, f_t, _) = sc.contact_forces(gap, normal, tang_disp);
        for (i, &v) in f_n.iter().enumerate() {
            assert_eq!(v, 0.0, "no normal force when separated at [{i}]");
        }
        let f_t_mag = (f_t[0] * f_t[0] + f_t[1] * f_t[1] + f_t[2] * f_t[2]).sqrt();
        assert!(f_t_mag < 1e-10, "no friction when separated");
    }
    #[test]
    fn test_sliding_contact_update_slip() {
        let mut sc = SlidingContact::new(1000.0, 500.0, 0.3, 2);
        sc.update_slip(0, [0.01, 0.02, 0.0]);
        assert!((sc.accumulated_slip[0][0] - 0.01).abs() < 1e-15);
        assert!((sc.accumulated_slip[0][1] - 0.02).abs() < 1e-15);
        sc.update_slip(0, [0.01, 0.0, 0.0]);
        assert!((sc.accumulated_slip[0][0] - 0.02).abs() < 1e-15);
        sc.reset_slip(0);
        assert_eq!(sc.accumulated_slip[0], [0.0; 3]);
    }
    #[test]
    fn test_aabb_from_points() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 0.0, 1.0], [2.0, 5.0, 2.0]];
        let aabb = Aabb::from_points(&points);
        assert_eq!(aabb.min, [1.0, 0.0, 1.0]);
        assert_eq!(aabb.max, [4.0, 5.0, 3.0]);
    }
    #[test]
    fn test_aabb_overlaps() {
        let a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [2.0, 2.0, 2.0],
        };
        let b = Aabb {
            min: [1.0, 1.0, 1.0],
            max: [3.0, 3.0, 3.0],
        };
        let c = Aabb {
            min: [5.0, 5.0, 5.0],
            max: [6.0, 6.0, 6.0],
        };
        assert!(a.overlaps(&b), "overlapping boxes");
        assert!(!a.overlaps(&c), "non-overlapping boxes");
    }
    #[test]
    fn test_aabb_contains_point() {
        let a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        assert!(a.contains_point([0.5, 0.5, 0.5]));
        assert!(!a.contains_point([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_aabb_volume() {
        let a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [2.0, 3.0, 4.0],
        };
        assert!((a.volume() - 24.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_expand() {
        let mut a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        a.expand(0.5);
        assert_eq!(a.min, [-0.5, -0.5, -0.5]);
        assert_eq!(a.max, [1.5, 1.5, 1.5]);
    }
    #[test]
    fn test_aabb_broad_phase() {
        let aabbs_a = vec![
            Aabb {
                min: [0.0, 0.0, 0.0],
                max: [1.0, 1.0, 1.0],
            },
            Aabb {
                min: [5.0, 5.0, 5.0],
                max: [6.0, 6.0, 6.0],
            },
        ];
        let aabbs_b = vec![Aabb {
            min: [0.5, 0.5, 0.5],
            max: [1.5, 1.5, 1.5],
        }];
        let pairs = AcceleratedContactSearch::aabb_broad_phase(&aabbs_a, &aabbs_b);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 0));
    }
    #[test]
    fn test_bucket_search() {
        let nodes_a: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let nodes_b: &[[f64; 3]] = &[[0.1, 0.1, 0.1], [5.0, 5.0, 5.0]];
        let pairs = AcceleratedContactSearch::bucket_search(nodes_a, nodes_b, 1.0, 0.5);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 0));
    }
    #[test]
    fn test_bucket_search_matches_brute_force() {
        let nodes_a: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let nodes_b: &[[f64; 3]] = &[[0.3, 0.0, 0.0], [1.2, 0.0, 0.0]];
        let cutoff = 0.5;
        let brute = ContactSearch::brute_force_search(nodes_a, nodes_b, cutoff);
        let bucket = AcceleratedContactSearch::bucket_search(nodes_a, nodes_b, 1.0, cutoff);
        assert_eq!(
            brute.candidate_pairs.len(),
            bucket.len(),
            "brute={} bucket={}",
            brute.candidate_pairs.len(),
            bucket.len()
        );
    }
    #[test]
    fn test_large_deformation_gap_updated() {
        let node = [0.0, 0.0, -1.5];
        let normal = [0.0, 0.0, 1.0];
        let gap = gap_function_node_to_plane(node, normal, 0.0);
        assert!(gap < 0.0, "node below plane: gap should be negative: {gap}");
    }
    #[test]
    fn test_large_deformation_contact_force() {
        let pc = PenaltyContact {
            stiffness: 500.0,
            friction_coeff: 0.2,
        };
        let gap = -0.5;
        let f = pc.normal_force(gap);
        assert!((f - 250.0).abs() < 1e-10, "f = {f}, expected 250");
    }
    #[test]
    fn test_self_contact_detection_nearby_nodes() {
        let nodes: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let pairs = self_contact_candidate_pairs(nodes, 0.2);
        assert_eq!(pairs.len(), 1, "one self-contact pair: {:?}", pairs);
        assert_eq!(pairs[0], (0, 1));
    }
    #[test]
    fn test_self_contact_no_self_pairs() {
        let nodes: &[[f64; 3]] = &[[0.0; 3], [1.0, 0.0, 0.0]];
        let pairs = self_contact_candidate_pairs(nodes, 0.5);
        assert!(pairs.is_empty(), "no self-contact pairs expected");
    }
    #[test]
    fn test_node_to_segment_at_midpoint() {
        let node = [1.0, 1.0, 0.0];
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [2.0, 0.0, 0.0];
        let (gap, normal) = gap_function_node_to_segment(node, p0, p1);
        assert!((gap - 1.0).abs() < 1e-12, "gap at midpoint: {gap}");
        assert!((normal[1] - 1.0).abs() < 1e-12, "normal y = {}", normal[1]);
    }
    #[test]
    fn test_node_to_segment_endpoint_clamping() {
        let node = [3.0, 1.0, 0.0];
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [2.0, 0.0, 0.0];
        let (gap, _) = gap_function_node_to_segment(node, p0, p1);
        assert!((gap - 2.0_f64.sqrt()).abs() < 1e-12, "gap = {gap}");
    }
    #[test]
    fn test_node_to_triangle_above_center() {
        let node = [0.5, 0.5, 1.0];
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let (gap, normal) = node_to_triangle_gap(node, v0, v1, v2);
        assert!((gap - 1.0).abs() < 1e-10, "gap = {gap}");
        assert!((normal[2] - 1.0).abs() < 1e-10, "normal z = {}", normal[2]);
    }
    #[test]
    fn test_node_to_triangle_below_plane() {
        let node = [0.25, 0.25, -0.5];
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let (gap, normal) = node_to_triangle_gap(node, v0, v1, v2);
        assert!(
            gap.abs() > 0.3,
            "below triangle: |gap| should be ~0.5, got {gap}"
        );
        assert!(
            normal[2].abs() > 0.9,
            "normal should be mostly z: n_z = {}",
            normal[2]
        );
    }
    #[test]
    fn test_node_to_triangle_on_plane() {
        let node = [0.25, 0.25, 0.0];
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let (gap, _normal) = node_to_triangle_gap(node, v0, v1, v2);
        assert!(gap.abs() < 1e-12, "on-plane gap should be zero: {gap}");
    }
    #[test]
    fn test_penalty_contact_stiffness_symmetry() {
        let pc = PenaltyContact {
            stiffness: 1000.0,
            friction_coeff: 0.3,
        };
        let normal = [0.0, 0.0, 1.0];
        let ke = penalty_contact_stiffness(pc.stiffness, normal);
        for (i, row) in ke.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert!(
                    (v - ke[j][i]).abs() < 1e-12,
                    "Stiffness not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_augmented_lagrangian_converges() {
        let mut alc = AugmentedLagrangianContact::new(100.0);
        for _ in 0..50 {
            alc.update_multipliers(-0.01, 0.0, 0.3);
        }
        assert!(
            alc.lambda_n < 0.0,
            "lambda_n should be negative: {}",
            alc.lambda_n
        );
    }
    #[test]
    fn test_contact_patch_area_positive() {
        let a = contact_patch_area(&[-0.01, -0.02, 0.005, 0.0], 1.0);
        assert!(a > 0.0, "contact patch area should be positive: {a}");
    }
    #[test]
    fn test_contact_patch_area_separated() {
        let a = contact_patch_area(&[0.01, 0.02, 0.03], 1.0);
        assert!(a.abs() < 1e-12, "no contact: area should be 0, got {a}");
    }
    #[test]
    fn test_hertz_contact_radius_positive() {
        let a = hertz_contact_radius(1000.0, 0.01, 200e9, 200e9, 0.3, 0.3);
        assert!(a > 0.0, "Hertz contact radius should be positive: {a}");
    }
    #[test]
    fn test_hertz_contact_pressure_positive() {
        let p0 = hertz_max_pressure(1000.0, 0.01, 200e9, 200e9, 0.3, 0.3);
        assert!(p0 > 0.0, "Hertz max pressure should be positive: {p0}");
    }
}
