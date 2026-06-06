//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BranchingCriterion, CrackTip3D, Fragment, PropagationCriterion, SplitResult};

/// Rankine (maximum principal stress) fracture criterion.
///
/// Fracture initiates when the maximum principal stress exceeds the tensile strength.
///
/// Input stress in Voigt notation: \[s11, s22, s33, s12, s23, s13\].
pub fn rankine_criterion(stress_voigt: [f64; 6], tensile_strength: f64) -> bool {
    let principals = principal_stresses_2d(stress_voigt);
    principals[0] > tensile_strength
}
/// Mohr-Coulomb fracture criterion.
///
/// `tau = c + sigma_n * tan(phi)`
///
/// Returns true if the criterion is exceeded (fracture initiates).
pub fn mohr_coulomb_criterion(
    stress_voigt: [f64; 6],
    cohesion: f64,
    friction_angle_rad: f64,
) -> bool {
    let principals = principal_stresses_2d(stress_voigt);
    let sigma1 = principals[0];
    let sigma3 = principals[2];
    let sin_phi = friction_angle_rad.sin();
    let cos_phi = friction_angle_rad.cos();
    let lhs = sigma1 * (1.0 + sin_phi) - sigma3 * (1.0 - sin_phi);
    let rhs = 2.0 * cohesion * cos_phi;
    lhs > rhs
}
/// Maximum shear stress (Tresca) criterion.
///
/// Fracture when max shear exceeds the shear strength.
pub fn tresca_criterion(stress_voigt: [f64; 6], shear_strength: f64) -> bool {
    let principals = principal_stresses_2d(stress_voigt);
    let tau_max = (principals[0] - principals[2]) / 2.0;
    tau_max > shear_strength
}
/// Von Mises (distortional energy) criterion.
///
/// Fracture when the von Mises stress exceeds the yield strength.
pub fn von_mises_criterion(stress_voigt: [f64; 6], yield_strength: f64) -> bool {
    let s = stress_voigt;
    let term1 = (s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2);
    let term2 = 6.0 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2));
    let vm = ((term1 + term2) / 2.0).sqrt();
    vm > yield_strength
}
/// Compute principal stresses from Voigt notation using the characteristic equation.
///
/// Returns sorted descending \[sigma1, sigma2, sigma3\].
pub(super) fn principal_stresses_2d(s: [f64; 6]) -> [f64; 3] {
    let sxx = s[0];
    let syy = s[1];
    let szz = s[2];
    let sxy = s[3];
    let syz = s[4];
    let sxz = s[5];
    let i1 = sxx + syy + szz;
    let i2 = sxx * syy + syy * szz + szz * sxx - sxy * sxy - syz * syz - sxz * sxz;
    let i3 = sxx * (syy * szz - syz * syz) - sxy * (sxy * szz - syz * sxz)
        + sxz * (sxy * syz - syy * sxz);
    let p3 = (i1 * i1 - 3.0 * i2) / 9.0;
    let q3 = (2.0 * i1.powi(3) - 9.0 * i1 * i2 + 27.0 * i3) / 54.0;
    let p3_safe = p3.max(0.0);
    let p_sqrt = p3_safe.sqrt();
    let val = if p3_safe > 1e-30 {
        (q3 / (p3_safe * p_sqrt)).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let theta = val.acos();
    let shift = i1 / 3.0;
    let r = 2.0 * p_sqrt;
    let mut eigs = [
        shift + r * (theta / 3.0).cos(),
        shift + r * ((theta + 2.0 * std::f64::consts::PI) / 3.0).cos(),
        shift + r * ((theta + 4.0 * std::f64::consts::PI) / 3.0).cos(),
    ];
    eigs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}
/// Griffith energy release rate for a crack of half-length `a` in a plate.
///
/// `G = pi * sigma^2 * a / E`
///
/// where `sigma` is the applied stress and `E` is Young's modulus.
/// For plane strain, multiply by `(1 - nu^2)`.
pub fn griffith_energy_release_rate(
    stress: f64,
    crack_half_length: f64,
    youngs_modulus: f64,
    plane_strain: bool,
    poisson_ratio: f64,
) -> f64 {
    if youngs_modulus.abs() < 1e-30 {
        return 0.0;
    }
    let g = std::f64::consts::PI * stress * stress * crack_half_length / youngs_modulus;
    if plane_strain {
        g * (1.0 - poisson_ratio * poisson_ratio)
    } else {
        g
    }
}
/// Stress intensity factor (mode I) for a central crack in infinite plate.
///
/// `K_I = sigma * sqrt(pi * a)`
pub fn stress_intensity_factor_mode1(stress: f64, crack_half_length: f64) -> f64 {
    stress * (std::f64::consts::PI * crack_half_length).sqrt()
}
/// Check if crack propagation occurs based on the Griffith criterion.
///
/// Crack grows when `G >= G_c` (critical energy release rate).
pub fn griffith_propagation_check(
    stress: f64,
    crack_half_length: f64,
    youngs_modulus: f64,
    critical_energy_release: f64,
    plane_strain: bool,
    poisson_ratio: f64,
) -> bool {
    let g = griffith_energy_release_rate(
        stress,
        crack_half_length,
        youngs_modulus,
        plane_strain,
        poisson_ratio,
    );
    g >= critical_energy_release
}
/// Split a node in a mesh to simulate crack opening.
///
/// Creates a duplicate of `node_idx` and reassigns edges on one side of the
/// crack plane to point to the new node.
///
/// `crack_normal` defines the crack plane normal at the node.
/// Edges whose opposite node lies on the positive side of the plane use the new node.
///
/// Returns a `SplitResult` with the new node index and reassigned edge indices.
pub fn split_node(
    positions: &mut Vec<[f64; 3]>,
    velocities: &mut Vec<[f64; 3]>,
    masses: &mut Vec<f64>,
    edges: &mut [(usize, usize)],
    node_idx: usize,
    crack_normal: [f64; 3],
) -> SplitResult {
    let pos = positions[node_idx];
    let vel = velocities[node_idx];
    let mass = masses[node_idx];
    let new_idx = positions.len();
    positions.push(pos);
    velocities.push(vel);
    masses.push(mass);
    let mut reassigned = Vec::new();
    for (edge_idx, edge) in edges.iter_mut().enumerate() {
        let other = if edge.0 == node_idx {
            edge.1
        } else if edge.1 == node_idx {
            edge.0
        } else {
            continue;
        };
        let dx = positions[other][0] - pos[0];
        let dy = positions[other][1] - pos[1];
        let dz = positions[other][2] - pos[2];
        let side = dx * crack_normal[0] + dy * crack_normal[1] + dz * crack_normal[2];
        if side > 0.0 {
            if edge.0 == node_idx {
                edge.0 = new_idx;
            } else {
                edge.1 = new_idx;
            }
            reassigned.push(edge_idx);
        }
    }
    SplitResult {
        original_node: node_idx,
        new_node: new_idx,
        reassigned_edges: reassigned,
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::fracture::CrackPath;

    use crate::fracture::DamageModel;
    use crate::fracture::FractureEdge;
    use crate::fracture::FractureMesh;
    use crate::fracture::FractureMode;
    use crate::fracture::FractureModel;

    use crate::fracture::griffith_energy_release_rate;

    use crate::fracture::stress_intensity_factor_mode1;

    #[test]
    fn test_fracture_edge_strain() {
        let edge = FractureEdge::new(0, 1, 1.0, 100.0);
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let strain = edge.current_strain(&positions);
        assert!(
            (strain - 0.5).abs() < 1e-10,
            "Expected strain 0.5, got {strain}"
        );
    }
    #[test]
    fn test_fracture_spring_force() {
        let stiffness = 200.0;
        let rest = 1.0;
        let edge = FractureEdge::new(0, 1, rest, stiffness);
        let pos_stretched = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let (fa, fb) = edge.spring_force(&pos_stretched);
        assert!(fa[0] > 0.0);
        assert!(fb[0] < 0.0);
        let expected_mag = stiffness * 0.5;
        assert!((fa[0] - expected_mag).abs() < 1e-8);
        let pos_compressed = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let (fa2, fb2) = edge.spring_force(&pos_compressed);
        assert!(fa2[0] < 0.0);
        assert!(fb2[0] > 0.0);
    }
    #[test]
    fn test_fracture_breaks_at_threshold() {
        let model = FractureModel {
            strain_threshold: 0.3,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut edge = FractureEdge::new(0, 1, 1.0, 100.0);
        let pos_safe = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]];
        let broke = edge.check_fracture(&pos_safe, &model);
        assert!(!broke);
        assert!(!edge.broken);
        let pos_broken = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let broke2 = edge.check_fracture(&pos_broken, &model);
        assert!(broke2);
        assert!(edge.broken);
        let (fa, fb) = edge.spring_force(&pos_broken);
        assert_eq!(fa, [0.0; 3]);
        assert_eq!(fb, [0.0; 3]);
    }
    #[test]
    fn test_fracture_mesh_gravity() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        let top = mesh.add_node([0.0, 1.0, 0.0], 1.0);
        let bot = mesh.add_node([0.0, 0.0, 0.0], 1.0);
        mesh.add_edge(top, bot, 500.0);
        mesh.pin_node(top);
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 60.0;
        let initial_y = mesh.positions[bot][1];
        for _ in 0..60 {
            mesh.step(dt, gravity);
        }
        let final_y = mesh.positions[bot][1];
        assert!(final_y < initial_y);
    }
    #[test]
    fn test_fracture_components_before() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        let a = mesh.add_node([0.0, 0.0, 0.0], 1.0);
        let b = mesh.add_node([1.0, 0.0, 0.0], 1.0);
        let c = mesh.add_node([2.0, 0.0, 0.0], 1.0);
        mesh.add_edge(a, b, 100.0);
        mesh.add_edge(b, c, 100.0);
        let components = mesh.connected_components();
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].len(), 3);
    }
    #[test]
    fn test_fracture_components_after() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        let a = mesh.add_node([0.0, 0.0, 0.0], 1.0);
        let b = mesh.add_node([1.0, 0.0, 0.0], 1.0);
        let c = mesh.add_node([2.0, 0.0, 0.0], 1.0);
        mesh.add_edge(a, b, 100.0);
        let mid = mesh.add_edge(b, c, 100.0);
        mesh.edges[mid].broken = true;
        let components = mesh.connected_components();
        assert_eq!(components.len(), 2);
    }
    #[test]
    fn test_fracture_step_no_panic() {
        let model = FractureModel {
            strain_threshold: 0.5,
            stress_threshold: 1000.0,
            mode: FractureMode::Combined,
        };
        let mut mesh = FractureMesh::new(model);
        let mut node_idx = [[0usize; 3]; 3];
        for (row, node_row) in node_idx.iter_mut().enumerate() {
            for (col, slot) in node_row.iter_mut().enumerate() {
                let x = col as f64;
                let y = row as f64;
                *slot = mesh.add_node([x, y, 0.0], 1.0);
            }
        }
        for &nid in node_idx[2].iter() {
            mesh.pin_node(nid);
        }
        for row_arr in node_idx.iter() {
            for col in 0..2 {
                mesh.add_edge(row_arr[col], row_arr[col + 1], 200.0);
            }
        }
        for (row, row_arr) in node_idx[..2].iter().enumerate() {
            for (col, &nid) in row_arr.iter().enumerate() {
                mesh.add_edge(nid, node_idx[row + 1][col], 200.0);
            }
        }
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 120.0;
        for _ in 0..100 {
            mesh.step(dt, gravity);
        }
        let _ = mesh.connected_components();
    }
    #[test]
    fn test_rankine_criterion() {
        let stress = [500.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(
            rankine_criterion(stress, 400.0),
            "Should fracture at 500 > 400"
        );
        assert!(
            !rankine_criterion(stress, 600.0),
            "Should not fracture at 500 < 600"
        );
    }
    #[test]
    fn test_mohr_coulomb() {
        let stress = [300.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let phi = 30.0_f64.to_radians();
        assert!(mohr_coulomb_criterion(stress, 10.0, phi));
        assert!(!mohr_coulomb_criterion(stress, 10000.0, phi));
    }
    #[test]
    fn test_tresca_criterion() {
        let stress = [200.0, 0.0, -100.0, 0.0, 0.0, 0.0];
        assert!(tresca_criterion(stress, 100.0), "tau_max=150 > 100");
        assert!(!tresca_criterion(stress, 200.0), "tau_max=150 < 200");
    }
    #[test]
    fn test_von_mises_criterion() {
        let stress = [300.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(von_mises_criterion(stress, 200.0), "VM should exceed 200");
        assert!(
            !von_mises_criterion(stress, 400.0),
            "VM should not exceed 400"
        );
    }
    #[test]
    fn test_griffith_energy_release() {
        let g = griffith_energy_release_rate(100.0, 0.01, 200e9, false, 0.3);
        assert!(g > 0.0, "Energy release rate should be positive");
        let g_ps = griffith_energy_release_rate(100.0, 0.01, 200e9, true, 0.3);
        assert!(g_ps > g * 0.9, "Plane strain G should be similar");
        assert!((g_ps - g * (1.0 - 0.3 * 0.3)).abs() < 1e-15);
    }
    #[test]
    fn test_stress_intensity_factor() {
        let k = stress_intensity_factor_mode1(100.0, 0.01);
        let expected = 100.0 * (std::f64::consts::PI * 0.01).sqrt();
        assert!(
            (k - expected).abs() < 1e-10,
            "K_I mismatch: {k} vs {expected}"
        );
    }
    #[test]
    fn test_damage_accumulation() {
        let mut dm = DamageModel::new(0.5, 0.1, 1000.0);
        assert!(!dm.is_broken());
        assert!((dm.effective_stiffness() - 1000.0).abs() < 1e-10);
        dm.accumulate(0.05);
        assert!((dm.damage - 0.0).abs() < 1e-10);
        dm.accumulate(0.3);
        assert!((dm.damage - 0.1).abs() < 1e-10);
        assert!((dm.effective_stiffness() - 900.0).abs() < 1e-10);
        dm.accumulate(10.0);
        assert!(dm.is_broken());
        assert!((dm.effective_stiffness() - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_crack_path() {
        let mut crack = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(crack.segment_count(), 0);
        assert!((crack.total_length() - 0.0).abs() < 1e-10);
        crack.propagate(1.0);
        assert_eq!(crack.segment_count(), 1);
        assert!((crack.total_length() - 1.0).abs() < 1e-10);
        let tip = crack.tip();
        assert!((tip[0] - 1.0).abs() < 1e-10);
        crack.propagate(0.5);
        assert_eq!(crack.segment_count(), 2);
        assert!((crack.total_length() - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_crack_direction_change() {
        let mut crack = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        crack.propagate(1.0);
        crack.propagate_with_direction(1.0, [0.0, 1.0, 0.0]);
        let tip = crack.tip();
        assert!((tip[0] - 1.0).abs() < 1e-10, "X should remain 1.0");
        assert!((tip[1] - 1.0).abs() < 1e-10, "Y should be 1.0");
    }
    #[test]
    fn test_split_node() {
        let mut positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let mut velocities = vec![[0.0; 3]; 3];
        let mut masses = vec![1.0; 3];
        let mut edges = vec![(0, 1), (0, 2)];
        let result = split_node(
            &mut positions,
            &mut velocities,
            &mut masses,
            &mut edges,
            0,
            [1.0, 0.0, 0.0],
        );
        assert_eq!(result.new_node, 3);
        assert_eq!(positions.len(), 4);
        assert_eq!(result.reassigned_edges.len(), 1);
        assert!(edges[1] == (0, 2) || edges[1] == (2, 0));
    }
    #[test]
    fn test_elastic_energy() {
        let edge = FractureEdge::new(0, 1, 1.0, 200.0);
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let e = edge.elastic_energy(&positions);
        assert!(
            (e - 25.0).abs() < 1e-10,
            "Elastic energy should be 25, got {e}"
        );
    }
    #[test]
    fn test_edge_damage_accumulation() {
        let mut edge = FractureEdge::new(0, 1, 1.0, 100.0);
        let positions = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]];
        let broke = edge.accumulate_damage(&positions, 1.0, 0.1);
        assert!(!broke);
        assert!((edge.damage - 0.1).abs() < 1e-10);
        assert!((edge.effective_stiffness() - 90.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_total_energy() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        let a = mesh.add_node([0.0, 0.0, 0.0], 1.0);
        let b = mesh.add_node([1.5, 0.0, 0.0], 1.0);
        mesh.add_edge(a, b, 100.0);
        let e = mesh.total_elastic_energy();
        assert!(e.abs() < 1e-10, "Energy at rest should be 0, got {e}");
    }
    #[test]
    fn test_max_strain() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        let a = mesh.add_node([0.0, 0.0, 0.0], 1.0);
        let b = mesh.add_node([1.0, 0.0, 0.0], 1.0);
        mesh.add_edge(a, b, 100.0);
        assert!(mesh.max_strain().abs() < 1e-10);
    }
    #[test]
    fn test_damped_step() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        let top = mesh.add_node([0.0, 1.0, 0.0], 1.0);
        let bot = mesh.add_node([0.0, 0.0, 0.0], 1.0);
        mesh.add_edge(top, bot, 500.0);
        mesh.pin_node(top);
        for _ in 0..100 {
            mesh.step_damped(1.0 / 60.0, [0.0, -9.81, 0.0], 0.01);
        }
        assert!(mesh.positions[bot][1] < 0.0);
    }
    #[test]
    fn test_griffith_propagation() {
        assert!(griffith_propagation_check(
            1000.0, 0.1, 200e9, 1e-6, false, 0.3
        ));
        assert!(!griffith_propagation_check(
            1.0, 0.0001, 200e9, 100.0, false, 0.3
        ));
    }
    #[test]
    fn test_kinetic_energy() {
        let model = FractureModel {
            strain_threshold: 10.0,
            stress_threshold: f64::MAX,
            mode: FractureMode::StrainBased,
        };
        let mut mesh = FractureMesh::new(model);
        mesh.add_node([0.0, 0.0, 0.0], 2.0);
        mesh.velocities[0] = [3.0, 4.0, 0.0];
        let ke = mesh.total_kinetic_energy();
        assert!((ke - 25.0).abs() < 1e-10, "KE should be 25, got {ke}");
    }
}
/// Compute the Mode-I stress intensity factor at a crack tip using the
/// finite-body correction formula (Irwin/Tada):
///
/// `K_I = σ √(π a) F(a/W)`
///
/// where `F(a/W)` is the geometry correction factor (edge-cracked plate):
/// `F = 1.12 - 0.23(a/W) + 10.56(a/W)^2 - 21.74(a/W)^3 + 30.42(a/W)^4`
///
/// # Arguments
/// * `stress` – remote applied stress (Pa)
/// * `crack_half_length` – half crack length a (m)
/// * `specimen_width` – specimen width W (m)
pub fn stress_intensity_3d(stress: f64, crack_half_length: f64, specimen_width: f64) -> f64 {
    let a = crack_half_length;
    let ratio = (a / specimen_width.max(f64::EPSILON)).min(1.0);
    let f =
        1.12 - 0.23 * ratio + 10.56 * ratio.powi(2) - 21.74 * ratio.powi(3) + 30.42 * ratio.powi(4);
    stress * (std::f64::consts::PI * a).sqrt() * f
}
/// Determine whether the crack should branch given dynamic stress intensity
/// factors and crack speed.
///
/// Returns the two branch directions if bifurcation occurs, otherwise `None`.
///
/// The branching angle follows the maximum circumferential stress criterion:
/// θ_branch ≈ ±arccos(3 K_II² / (K_I² + K_II²)) simplified.
pub fn crack_branching_3d(
    k1: f64,
    k2: f64,
    crack_speed: f64,
    wave_speed: f64,
) -> Option<([f64; 3], [f64; 3])> {
    let speed_frac = crack_speed / wave_speed.max(f64::EPSILON);
    if speed_frac < 0.35 {
        return None;
    }
    let mix = (k2 / k1.abs().max(f64::EPSILON)).abs();
    if mix < 0.2 {
        return None;
    }
    let theta = (k2 / k1.abs().max(f64::EPSILON)).atan();
    let branch_angle = theta.abs() + std::f64::consts::PI / 8.0;
    let branch1 = [branch_angle.cos(), branch_angle.sin(), 0.0];
    let branch2 = [branch_angle.cos(), -branch_angle.sin(), 0.0];
    Some((branch1, branch2))
}
/// Count connected fragments in a fractured mesh using union-find.
///
/// Each node starts in its own component; intact (unbroken) edges merge
/// their endpoints.  Returns the number of distinct connected components.
pub fn count_fragments(n_nodes: usize, edges: &[(usize, usize)], broken: &[bool]) -> usize {
    let mut parent: Vec<usize> = (0..n_nodes).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    for (idx, &(a, b)) in edges.iter().enumerate() {
        if idx < broken.len() && broken[idx] {
            continue;
        }
        let ra = find(&mut parent, a);
        let rb = find(&mut parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }
    let mut roots: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for i in 0..n_nodes {
        roots.insert(find(&mut parent, i));
    }
    roots.len()
}
/// Check whether a crack should arrest using the Ritchie-Knott-Rice (RKR)
/// crack arrest criterion.
///
/// Arrest occurs when the stress intensity factor drops below the arrest
/// toughness `K_Ia`.
///
/// Returns `true` if the crack arrests.
pub fn crack_arrest_rkr(k_current: f64, k_arrest_toughness: f64) -> bool {
    k_current < k_arrest_toughness
}
/// Crack arrest based on a thermal blunting model.
///
/// At high temperatures the material toughness increases.  Returns true if the
/// effective toughness at temperature `T` (°C) exceeds the current SIF.
///
/// Uses a simple Arrhenius-like toughness: K_c(T) = K_c0 * exp(α * T / T_ref).
pub fn crack_arrest_thermal(
    k_current: f64,
    k_c0: f64,
    alpha: f64,
    temperature: f64,
    t_ref: f64,
) -> bool {
    let k_c_eff = k_c0 * (alpha * temperature / t_ref.max(f64::EPSILON)).exp();
    k_current < k_c_eff
}
/// Compute the limiting dynamic crack speed (Freund 1990).
///
/// The dynamic SIF is reduced by a factor (1 - v/c_R) where v is crack speed
/// and c_R is the Rayleigh wave speed.
///
/// `K_dyn = K_static * (1 - v / c_R)`
///
/// The crack cannot propagate faster than `c_R`.
pub fn dynamic_sif(k_static: f64, crack_speed: f64, rayleigh_speed: f64) -> f64 {
    let v = crack_speed.min(rayleigh_speed).max(0.0);
    k_static * (1.0 - v / rayleigh_speed.max(f64::EPSILON))
}
/// Estimate the terminal crack speed given a far-field SIF `k_inf` and the
/// crack-initiation toughness `k_0`:
///
/// `v_terminal ≈ c_R * (1 - (K_0 / K_inf)²)` (Broberg 1999 approximation).
pub fn terminal_crack_speed(k_inf: f64, k_0: f64, rayleigh_speed: f64) -> f64 {
    if k_inf.abs() < f64::EPSILON {
        return 0.0;
    }
    let ratio = k_0 / k_inf;
    rayleigh_speed * (1.0 - ratio * ratio).max(0.0)
}
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
#[cfg(test)]
mod tests_fracture_ext {

    use crate::fracture::CrackTip3D;

    use crate::fracture::count_fragments;
    use crate::fracture::crack_arrest_rkr;
    use crate::fracture::crack_arrest_thermal;

    use crate::fracture::crack_branching_3d;

    use crate::fracture::dynamic_sif;

    use crate::fracture::stress_intensity_3d;

    use crate::fracture::terminal_crack_speed;
    #[test]
    fn test_crack_tip_3d_advance() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        tip.advance(2.5);
        assert!((tip.position[0] - 2.5).abs() < 1e-10);
        assert!((tip.length - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_crack_tip_3d_arrested() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        tip.arrested = true;
        tip.advance(5.0);
        assert!(
            tip.position[0].abs() < 1e-12,
            "arrested tip should not move"
        );
    }
    #[test]
    fn test_crack_tip_3d_steer() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        tip.steer([0.0, 1.0, 0.0]);
        assert!((tip.direction[1] - 1.0).abs() < 1e-10);
        assert!(tip.direction[0].abs() < 1e-10);
    }
    #[test]
    fn test_stress_intensity_3d_pure_tension() {
        let k = stress_intensity_3d(100e6, 0.01, 0.1);
        let approx = 100e6_f64 * (std::f64::consts::PI * 0.01_f64).sqrt() * 1.12;
        assert!(
            (k - approx).abs() / approx < 0.30,
            "K_I = {k}, approx = {approx}"
        );
    }
    #[test]
    fn test_crack_branching_3d_slow_no_branch() {
        let result = crack_branching_3d(1e6, 3e5, 100.0, 3000.0);
        assert!(result.is_none(), "slow crack should not branch");
    }
    #[test]
    fn test_crack_branching_3d_fast_mixed_mode() {
        let result = crack_branching_3d(1e6, 4e5, 1200.0, 3000.0);
        assert!(result.is_some(), "fast mixed-mode crack should branch");
    }
    #[test]
    fn test_count_fragments_no_breaks() {
        let edges = vec![(0, 1), (1, 2), (2, 3)];
        let broken = vec![false, false, false];
        let n = count_fragments(4, &edges, &broken);
        assert_eq!(n, 1, "fully connected → 1 fragment");
    }
    #[test]
    fn test_count_fragments_all_broken() {
        let edges = vec![(0, 1), (1, 2), (2, 3)];
        let broken = vec![true, true, true];
        let n = count_fragments(4, &edges, &broken);
        assert_eq!(n, 4, "all broken → 4 fragments");
    }
    #[test]
    fn test_count_fragments_partial() {
        let edges = vec![(0, 1), (1, 2), (2, 3)];
        let broken = vec![false, true, false];
        let n = count_fragments(4, &edges, &broken);
        assert_eq!(n, 2);
    }
    #[test]
    fn test_crack_arrest_rkr_below_toughness() {
        assert!(crack_arrest_rkr(0.5e6, 1.0e6), "K < K_Ia → arrest");
        assert!(!crack_arrest_rkr(2.0e6, 1.0e6), "K > K_Ia → no arrest");
    }
    #[test]
    fn test_crack_arrest_thermal_high_temp() {
        assert!(crack_arrest_thermal(1e6, 1e6, 5.0, 500.0, 300.0));
    }
    #[test]
    fn test_dynamic_sif_zero_speed() {
        let k_dyn = dynamic_sif(1e6, 0.0, 3000.0);
        assert!((k_dyn - 1e6).abs() < 1.0, "zero speed → full static SIF");
    }
    #[test]
    fn test_dynamic_sif_at_rayleigh_speed() {
        let k_dyn = dynamic_sif(1e6, 3000.0, 3000.0);
        assert!(k_dyn.abs() < 1.0, "at Rayleigh speed → K_dyn ≈ 0");
    }
    #[test]
    fn test_terminal_crack_speed_high_k_inf() {
        let v = terminal_crack_speed(10e6, 1e6, 3000.0);
        assert!(v > 2900.0, "v_terminal should approach c_R: {v}");
    }
    #[test]
    fn test_terminal_crack_speed_zero_k_inf() {
        let v = terminal_crack_speed(0.0, 1e6, 3000.0);
        assert!(v.abs() < 1e-10, "zero K_inf → zero speed");
    }
}
/// Build a unit vector perpendicular to `n` using a stable Gram-Schmidt step.
pub(super) fn make_perpendicular(n: [f64; 3]) -> [f64; 3] {
    let ax: [f64; 3] = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let dot = ax[0] * n[0] + ax[1] * n[1] + ax[2] * n[2];
    let t = [ax[0] - dot * n[0], ax[1] - dot * n[1], ax[2] - dot * n[2]];
    let len = norm3(t);
    if len > 1e-12 {
        [t[0] / len, t[1] / len, t[2] / len]
    } else {
        [1.0, 0.0, 0.0]
    }
}
/// Propagate a single crack-front tip using the selected criterion.
///
/// # Arguments
/// * `tip`         – mutable crack-tip state (position, direction, length)
/// * `k1`          – Mode-I stress intensity factor (Pa√m)
/// * `k2`          – Mode-II stress intensity factor (Pa√m)
/// * `toughness`   – fracture toughness K_Ic (Pa√m)
/// * `da`          – crack advance increment (m)
/// * `criterion`   – which propagation rule to apply
///
/// The function advances `tip` by `da` in the direction determined by
/// the chosen criterion.  Returns `true` when the crack propagated, `false`
/// when it arrested (K_eq < toughness).
pub fn propagate_crack(
    tip: &mut CrackTip3D,
    k1: f64,
    k2: f64,
    toughness: f64,
    da: f64,
    criterion: PropagationCriterion,
) -> bool {
    if tip.arrested {
        return false;
    }
    let k_eq = match criterion {
        PropagationCriterion::MaxHoopStress => {
            let theta = hoop_stress_angle(k1, k2);
            let ch = (theta / 2.0).cos();
            let s = theta.sin();
            let base = ch * (k1 * ch * ch - 1.5 * k2 * s);
            base.abs()
        }
        PropagationCriterion::MinStrainEnergy => (k1 * k1 + k2 * k2).sqrt(),
        PropagationCriterion::CohesiveZone => k1.abs(),
    };
    if k_eq < toughness {
        tip.arrested = true;
        return false;
    }
    let new_dir = propagation_direction(k1, k2, tip.direction, criterion);
    tip.steer(new_dir);
    tip.advance(da);
    true
}
/// Compute the propagation direction vector for the given stress intensities
/// and current crack direction, according to `criterion`.
///
/// Returns a unit vector in 3-D.
pub fn propagation_direction(
    k1: f64,
    k2: f64,
    current_direction: [f64; 3],
    criterion: PropagationCriterion,
) -> [f64; 3] {
    let theta = match criterion {
        PropagationCriterion::MaxHoopStress => hoop_stress_angle(k1, k2),
        PropagationCriterion::MinStrainEnergy => min_strain_energy_angle(k1, k2),
        PropagationCriterion::CohesiveZone => 0.0,
    };
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let n = norm3(current_direction);
    let d = if n > 1e-12 {
        [
            current_direction[0] / n,
            current_direction[1] / n,
            current_direction[2] / n,
        ]
    } else {
        [1.0, 0.0, 0.0]
    };
    let t = make_perp_fracture(d);
    [
        cos_t * d[0] + sin_t * t[0],
        cos_t * d[1] + sin_t * t[1],
        cos_t * d[2] + sin_t * t[2],
    ]
}
/// Maximum hoop stress propagation angle θ_c (Williams 1957).
#[inline]
pub(super) fn hoop_stress_angle(k1: f64, k2: f64) -> f64 {
    if k2.abs() < 1e-30 {
        return 0.0;
    }
    let ratio = k1 / (4.0 * k2);
    let sign = if k2 >= 0.0 { 1.0_f64 } else { -1.0 };
    let inner = ratio * ratio + 0.5;
    2.0 * (ratio - sign * inner.sqrt()).atan()
}
/// Minimum strain energy density angle (Sih 1974).
#[inline]
pub(super) fn min_strain_energy_angle(k1: f64, k2: f64) -> f64 {
    if k1.abs() < 1e-30 {
        return 0.0;
    }
    let ratio = k2 / k1;
    let disc = (1.0 + 8.0 * ratio * ratio).sqrt();
    -2.0 * (ratio / (1.0 + disc)).atan()
}
/// Internal perpendicular helper (avoids name collision with `make_perpendicular`).
#[inline]
pub(super) fn make_perp_fracture(n: [f64; 3]) -> [f64; 3] {
    let ax = if n[0].abs() < 0.9 {
        [1.0_f64, 0.0, 0.0]
    } else {
        [0.0_f64, 1.0, 0.0]
    };
    let dot = ax[0] * n[0] + ax[1] * n[1] + ax[2] * n[2];
    let t = [ax[0] - dot * n[0], ax[1] - dot * n[1], ax[2] - dot * n[2]];
    let len = norm3(t);
    if len > 1e-12 {
        [t[0] / len, t[1] / len, t[2] / len]
    } else {
        [0.0, 1.0, 0.0]
    }
}
/// Free-function wrapper for crack branching.
///
/// Applies the default [`BranchingCriterion`] parameters.
///
/// Returns `None` if the crack does not branch; `Some((tip_a, tip_b))`
/// otherwise, where `tip_a` and `tip_b` are the two new crack tips.
pub fn crack_branching(
    tip: &CrackTip3D,
    k1: f64,
    k2: f64,
    k1c: f64,
) -> Option<(CrackTip3D, CrackTip3D)> {
    BranchingCriterion::default_params().crack_branching(tip, k1, k2, k1c)
}
/// Estimate the Rayleigh wave speed from material properties.
///
/// Approximation from Viktorov (1967):
/// `c_R ≈ c_S * (0.862 + 1.14 ν) / (1 + ν)`
///
/// where `c_S = sqrt(μ / ρ)` is the shear wave speed.
pub fn rayleigh_wave_speed(young: f64, poisson: f64, density: f64) -> f64 {
    if density < 1e-30 || young < 1e-30 {
        return 0.0;
    }
    let mu = young / (2.0 * (1.0 + poisson));
    let c_s = (mu / density).sqrt();
    let nu = poisson;
    c_s * (0.862 + 1.14 * nu) / (1.0 + nu)
}
/// Limiting crack speed as a fraction of the Rayleigh wave speed.
///
/// Empirical result (Stroh 1954, Freund 1990): cracks in brittle materials
/// seldom exceed `≈ 0.6 c_R` before branching.  Returns the theoretical
/// maximum in m/s.
pub fn limiting_crack_speed(young: f64, poisson: f64, density: f64) -> f64 {
    0.6 * rayleigh_wave_speed(young, poisson, density)
}
#[cfg(test)]
mod tests_crack_front {
    use super::*;
    use crate::fracture::BranchingCriterion;

    use crate::fracture::CrackFront;
    use crate::fracture::CrackFrontVertex;

    use crate::fracture::CrackTip3D;

    use crate::fracture::PropagationCriterion;

    use crate::fracture::crack_branching;

    use crate::fracture::griffith_energy_release_rate;

    use crate::fracture::limiting_crack_speed;
    use crate::fracture::make_perpendicular;

    use crate::fracture::norm3;
    use crate::fracture::propagation_direction;
    use crate::fracture::rayleigh_wave_speed;

    use crate::fracture::stress_intensity_factor_mode1;

    fn make_front(k1: f64, k2: f64) -> CrackFront {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        let v = CrackFrontVertex::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], k1, k2);
        front.add_vertex(v);
        front
    }
    #[test]
    fn test_crack_front_vertex_new_normalises_normal() {
        let v = CrackFrontVertex::new([0.0; 3], [3.0, 4.0, 0.0], 1e6, 0.0);
        let len = norm3(v.normal);
        assert!(
            (len - 1.0).abs() < 1e-12,
            "normal should be unit length, got {len}"
        );
    }
    #[test]
    fn test_max_hoop_stress_pure_mode1_is_zero() {
        let v = CrackFrontVertex::new([0.0; 3], [0.0, 1.0, 0.0], 2e6, 0.0);
        assert!(
            v.max_hoop_stress_angle().abs() < 1e-12,
            "pure mode-I → θ = 0"
        );
    }
    #[test]
    fn test_equivalent_sif_pure_mode1() {
        let v = CrackFrontVertex::new([0.0; 3], [0.0, 1.0, 0.0], 2e6, 0.0);
        assert!(
            (v.equivalent_sif() - 2e6).abs() < 1e-3,
            "K_eq should equal K_I for pure mode-I"
        );
    }
    #[test]
    fn test_crack_front_no_propagation_below_toughness() {
        let mut front = make_front(0.5e6, 0.0);
        let n = front.propagate();
        assert_eq!(n, 0, "sub-critical crack should not propagate");
    }
    #[test]
    fn test_crack_front_propagates_above_toughness() {
        let mut front = make_front(3e6, 0.0);
        let start = front.vertices[0].position;
        let n = front.propagate();
        assert_eq!(n, 1, "supercritical crack should propagate");
        let moved = (0..3).any(|d| (front.vertices[0].position[d] - start[d]).abs() > 1e-10);
        assert!(moved, "front vertex should have moved after propagation");
    }
    #[test]
    fn test_crack_front_propagation_increment_clamped() {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        let v = CrackFrontVertex::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1e9, 0.0);
        front.add_vertex(v);
        front.propagate();
        let pos = front.vertices[0].position;
        let moved = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        assert!(
            moved <= 1e-2 + 1e-12,
            "increment should be clamped to da_max: {moved}"
        );
    }
    #[test]
    fn test_crack_front_arrested_vertex_does_not_move() {
        let mut front = make_front(5e6, 0.0);
        front.vertices[0].arrested = true;
        let start = front.vertices[0].position;
        front.propagate();
        assert_eq!(
            front.vertices[0].position, start,
            "arrested vertex should not move"
        );
    }
    #[test]
    fn test_crack_front_branching_below_threshold_no_branches() {
        let front = make_front(1.5e6, 0.0);
        let branches = front.check_branching();
        assert!(branches.is_empty(), "below branch threshold → no branches");
    }
    #[test]
    fn test_crack_front_branching_above_threshold_two_branches() {
        let front = make_front(5e6, 0.0);
        let branches = front.check_branching();
        assert_eq!(branches.len(), 2, "should produce exactly 2 branches");
    }
    #[test]
    fn test_crack_front_branching_normals_are_unit() {
        let front = make_front(5e6, 1e6);
        let branches = front.check_branching();
        for (i, b) in branches.iter().enumerate() {
            let len = norm3(b.normal);
            assert!(
                (len - 1.0).abs() < 1e-10,
                "branch {i} normal should be unit: {len}"
            );
        }
    }
    #[test]
    fn test_apply_branches_increases_vertex_count() {
        let mut front = make_front(5e6, 0.0);
        let branches = front.check_branching();
        let n_before = front.vertices.len();
        front.apply_branches(branches);
        assert_eq!(front.vertices.len(), n_before + 2);
    }
    #[test]
    fn test_delete_fractured_elements_near_front() {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        let v = CrackFrontVertex::new([0.5, 0.5, 0.5], [0.0, 1.0, 0.0], 2e6, 0.0);
        front.add_vertex(v);
        let tet_indices: &[[usize; 4]] = &[[0, 1, 2, 3], [4, 5, 6, 7]];
        let positions: &[[f64; 3]] = &[
            [0.5, 0.5, 0.5],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [10.0, 10.0, 10.0],
            [11.0, 10.0, 10.0],
            [10.0, 11.0, 10.0],
            [10.0, 10.0, 11.0],
        ];
        let deleted = front.delete_fractured_elements(tet_indices, positions, 0.1);
        assert_eq!(deleted, 1, "one element should be deleted");
        assert!(front.is_element_deleted(0), "element 0 should be deleted");
        assert!(
            !front.is_element_deleted(1),
            "element 1 should not be deleted"
        );
    }
    #[test]
    fn test_delete_fractured_elements_no_deletion_when_far() {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        let v = CrackFrontVertex::new([100.0, 100.0, 100.0], [0.0, 1.0, 0.0], 2e6, 0.0);
        front.add_vertex(v);
        let tet_indices: &[[usize; 4]] = &[[0, 1, 2, 3]];
        let positions: &[[f64; 3]] = &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let deleted = front.delete_fractured_elements(tet_indices, positions, 0.5);
        assert_eq!(
            deleted, 0,
            "no element should be deleted when far from front"
        );
    }
    #[test]
    fn test_arrest_sub_critical_marks_low_k_vertices() {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        front.add_vertex(CrackFrontVertex::new([0.0; 3], [0.0, 1.0, 0.0], 0.3e6, 0.0));
        front.add_vertex(CrackFrontVertex::new(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            3e6,
            0.0,
        ));
        front.arrest_sub_critical();
        assert!(
            front.vertices[0].arrested,
            "sub-critical vertex should be arrested"
        );
        assert!(
            !front.vertices[1].arrested,
            "supercritical vertex should not be arrested"
        );
    }
    #[test]
    fn test_active_vertex_count() {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        front.add_vertex(CrackFrontVertex::new([0.0; 3], [0.0, 1.0, 0.0], 2e6, 0.0));
        front.add_vertex(CrackFrontVertex {
            position: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            k1: 2e6,
            k2: 0.0,
            arrested: true,
        });
        assert_eq!(front.active_vertex_count(), 1);
    }
    #[test]
    fn test_estimate_sif_nearest_node() {
        let mut front = CrackFront::new(1e6, 2e6, 1e-4, 1e-2);
        let v = CrackFrontVertex::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, 0.0);
        front.add_vertex(v);
        let positions = vec![[0.01_f64, 0.0, 0.0]];
        let stress_voigt = vec![[0.0_f64, 1e6, 0.0, 0.0, 0.0, 0.0]];
        front.estimate_sif_from_stress_field(&stress_voigt, &positions);
        assert!(
            front.vertices[0].k1 > 0.0,
            "K_I should be estimated > 0 for tensile stress"
        );
    }
    #[test]
    fn test_make_perpendicular_orthogonality() {
        for n in [
            [1.0_f64, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.577, 0.577, 0.577_f64],
        ] {
            let nl = norm3(n);
            let n_unit = [n[0] / nl, n[1] / nl, n[2] / nl];
            let t = make_perpendicular(n_unit);
            let dot = t[0] * n_unit[0] + t[1] * n_unit[1] + t[2] * n_unit[2];
            assert!(
                dot.abs() < 1e-10,
                "perpendicular dot product should be 0, got {dot}"
            );
            let tl = norm3(t);
            assert!(
                (tl - 1.0).abs() < 1e-10,
                "perpendicular should be unit length, got {tl}"
            );
        }
    }
    #[test]
    fn test_propagate_crack_max_hoop_stress_pure_mode1() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let k1 = 2e6;
        let k2 = 0.0;
        let toughness = 1e6;
        let da = 0.01;
        let propagated = propagate_crack(
            &mut tip,
            k1,
            k2,
            toughness,
            da,
            PropagationCriterion::MaxHoopStress,
        );
        assert!(propagated, "crack should propagate when K_I > toughness");
        assert!(
            (tip.position[0] - da).abs() < 1e-8,
            "crack should advance in X: pos={:?}",
            tip.position
        );
        assert!((tip.length - da).abs() < 1e-10);
    }
    #[test]
    fn test_propagate_crack_arrests_below_toughness() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let propagated = propagate_crack(
            &mut tip,
            0.5e6,
            0.0,
            1e6,
            0.01,
            PropagationCriterion::MaxHoopStress,
        );
        assert!(!propagated, "crack should arrest when K_eq < toughness");
        assert!(tip.arrested, "tip should be arrested");
        assert!(
            tip.position[0].abs() < 1e-15,
            "arrested crack should not move"
        );
    }
    #[test]
    fn test_propagate_crack_min_strain_energy() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let propagated = propagate_crack(
            &mut tip,
            2e6,
            5e5,
            1e6,
            0.01,
            PropagationCriterion::MinStrainEnergy,
        );
        assert!(
            propagated,
            "crack should propagate (min-strain-energy criterion)"
        );
        assert!(tip.length > 0.0);
    }
    #[test]
    fn test_propagate_crack_cohesive_zone() {
        let mut tip = CrackTip3D::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let propagated = propagate_crack(
            &mut tip,
            3e6,
            0.0,
            1e6,
            0.005,
            PropagationCriterion::CohesiveZone,
        );
        assert!(propagated, "cohesive zone criterion: should propagate");
        assert!((tip.length - 0.005).abs() < 1e-10);
    }
    #[test]
    fn test_propagation_direction_pure_mode1_straight() {
        let dir = propagation_direction(
            1e6,
            0.0,
            [1.0, 0.0, 0.0],
            PropagationCriterion::MaxHoopStress,
        );
        assert!(
            (dir[0] - 1.0).abs() < 1e-6,
            "pure Mode-I should propagate straight: dir={dir:?}"
        );
    }
    #[test]
    fn test_propagation_direction_unit_length() {
        let dir = propagation_direction(
            1e6,
            3e5,
            [1.0, 0.0, 0.0],
            PropagationCriterion::MaxHoopStress,
        );
        let len = norm3(dir);
        assert!(
            (len - 1.0).abs() < 1e-10,
            "propagation direction should be unit vector: {len}"
        );
    }
    #[test]
    fn test_crack_branching_below_k_threshold_no_branch() {
        let tip = CrackTip3D::new([0.0; 3], [1.0, 0.0, 0.0]);
        let result = crack_branching(&tip, 1e6, 0.0, 1e6);
        assert!(result.is_none(), "K_eq < 2 * K_Ic should not branch");
    }
    #[test]
    fn test_crack_branching_high_k_mixed_mode_branches() {
        let tip = CrackTip3D::new([0.0; 3], [1.0, 0.0, 0.0]);
        let result = crack_branching(&tip, 3e6, 2e6, 1e6);
        assert!(result.is_some(), "large mixed-mode K should branch");
        let (a, b) = result.unwrap();
        let same = (0..3).all(|d| (a.direction[d] - b.direction[d]).abs() < 1e-10);
        assert!(!same, "two branches should have different directions");
    }
    #[test]
    fn test_crack_branching_branch_tips_at_same_origin() {
        let tip = CrackTip3D::new([1.0, 2.0, 3.0], [1.0, 0.0, 0.0]);
        if let Some((a, b)) = crack_branching(&tip, 3e6, 2e6, 1e6) {
            for d in 0..3 {
                assert!((a.position[d] - tip.position[d]).abs() < 1e-12);
                assert!((b.position[d] - tip.position[d]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_branching_criterion_custom_params() {
        let criterion = BranchingCriterion {
            angle_threshold: 0.0,
            k_ratio_threshold: 1.0,
            branch_half_angle: std::f64::consts::PI / 6.0,
        };
        let tip = CrackTip3D::new([0.0; 3], [1.0, 0.0, 0.0]);
        let result = criterion.crack_branching(&tip, 2e6, 1e6, 1e6);
        assert!(result.is_some(), "relaxed criterion should allow branching");
    }
    #[test]
    fn test_rayleigh_wave_speed_steel() {
        let c_r = rayleigh_wave_speed(200e9, 0.3, 7800.0);
        assert!(
            c_r > 2500.0 && c_r < 3500.0,
            "Rayleigh speed for steel should be ~2900 m/s, got {c_r}"
        );
    }
    #[test]
    fn test_rayleigh_wave_speed_zero_density() {
        let c_r = rayleigh_wave_speed(200e9, 0.3, 0.0);
        assert_eq!(c_r, 0.0, "zero density → zero speed");
    }
    #[test]
    fn test_limiting_crack_speed_less_than_rayleigh() {
        let c_r = rayleigh_wave_speed(200e9, 0.3, 7800.0);
        let c_lim = limiting_crack_speed(200e9, 0.3, 7800.0);
        assert!(
            c_lim < c_r,
            "limiting crack speed should be less than Rayleigh speed"
        );
        assert!((c_lim - 0.6 * c_r).abs() < 1e-6);
    }
    #[test]
    fn test_energy_release_rate_consistent_with_propagation() {
        let stress = 500.0_f64;
        let a = 0.05_f64;
        let e = 200e9_f64;
        let g = griffith_energy_release_rate(stress, a, e, false, 0.3);
        assert!(g > 0.0, "energy release rate should be positive");
        let k1 = stress_intensity_factor_mode1(stress, a);
        let g_from_k = k1 * k1 / e;
        let rel = (g - g_from_k).abs() / (g + 1e-30);
        assert!(
            rel < 0.01,
            "G from Griffith ≈ K_I^2/E: G={g}, K²/E={g_from_k}, rel={rel}"
        );
    }
}
/// Estimate the number of fragments from a dynamic impact using the
/// Grady–Kipp energy model.
///
/// `N ≈ (ρ c² ε̇²) / (2 K_Ic²)`  (approximate, dimensionless count per volume).
///
/// In practice the number of macro fragments scales with impact energy / G_c:
///
/// `N_frags = max(1, round(E_impact / G_c))`
///
/// where `G_c = K_Ic² / E_young` (plane stress) is the fracture energy per
/// unit volume times a characteristic length.
pub fn dynamic_fragment_count(
    impact_energy: f64,
    k1c: f64,
    young: f64,
    characteristic_length: f64,
) -> usize {
    if young < 1e-30 || k1c < 1e-30 {
        return 1;
    }
    let g_c = k1c * k1c / young;
    let energy_per_fragment = g_c * characteristic_length * characteristic_length;
    if energy_per_fragment < 1e-30 {
        return 1;
    }
    let n = (impact_energy / energy_per_fragment).round() as usize;
    n.max(1)
}
/// Estimate the average fragment size (m) using the Mott fragmentation formula.
///
/// For a ring or rod under tension expanding at strain rate `strain_rate` (1/s):
///
/// `L_frag ≈ K * (K_Ic / (ρ c²))^(2/3) * (C / strain_rate)^(2/3)`
///
/// Simplified form used here:
///
/// `L_frag = (2 K_Ic²) / (ρ c² strain_rate²)`^(1/3)
///
/// where `c = sqrt(E/ρ)` is the longitudinal wave speed.
pub fn mott_fragment_size(k1c: f64, density: f64, young: f64, strain_rate: f64) -> f64 {
    if density < 1e-30 || young < 1e-30 || strain_rate.abs() < 1e-30 {
        return f64::MAX;
    }
    let c2 = young / density;
    let numerator = 2.0 * k1c * k1c;
    let denominator = density * c2 * strain_rate * strain_rate;
    if denominator < 1e-60 {
        return f64::MAX;
    }
    (numerator / denominator).powf(1.0 / 3.0)
}
/// Extract individual fragments from a fractured mesh.
///
/// Uses BFS over non-broken edges to find connected components, then builds
/// one [`Fragment`] per component.
pub fn extract_fragments(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    masses: &[f64],
    edges: &[(usize, usize)],
    broken: &[bool],
) -> Vec<Fragment> {
    let n = positions.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (idx, &(a, b)) in edges.iter().enumerate() {
        if idx < broken.len() && broken[idx] {
            continue;
        }
        adj[a].push(b);
        adj[b].push(a);
    }
    let mut visited = vec![false; n];
    let mut fragments = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(node) = queue.pop_front() {
            component.push(node);
            for &nb in &adj[node] {
                if !visited[nb] {
                    visited[nb] = true;
                    queue.push_back(nb);
                }
            }
        }
        fragments.push(Fragment::from_nodes(
            component, positions, velocities, masses,
        ));
    }
    fragments
}
/// Stress intensity factor for an edge crack in a finite-width plate.
///
/// Uses the Tada (1985) correction:
/// `K_I = σ √(π a) * F(a/W)`
///
/// where `F(a/W) ≈ 1.12 - 0.231 (a/W) + 10.55 (a/W)^2 - 21.72 (a/W)^3 + 30.39 (a/W)^4`
pub fn sif_edge_crack_finite_plate(sigma: f64, crack_length: f64, plate_width: f64) -> f64 {
    let ratio = (crack_length / plate_width.max(1e-30)).min(1.0);
    let f = 1.12 - 0.231 * ratio + 10.55 * ratio.powi(2) - 21.72 * ratio.powi(3)
        + 30.39 * ratio.powi(4);
    sigma * (std::f64::consts::PI * crack_length).sqrt() * f
}
/// Griffith critical stress for crack initiation.
///
/// `σ_c = √(2 E G_c / (π a))` (plane stress)
///
/// Returns the minimum far-field stress required to propagate a crack of
/// half-length `a`.
pub fn griffith_critical_stress(young: f64, fracture_energy: f64, crack_half_length: f64) -> f64 {
    if crack_half_length < 1e-30 {
        return f64::MAX;
    }
    (2.0 * young * fracture_energy / (std::f64::consts::PI * crack_half_length)).sqrt()
}
/// Crack opening displacement (COD) at the crack tip for mode-I loading.
///
/// `δ_t = 4 K_I² / (π E σ_ys)` (Dugdale strip-yield model)
///
/// where `σ_ys` is the yield stress.
pub fn crack_opening_displacement(k1: f64, young: f64, yield_stress: f64) -> f64 {
    if young < 1e-30 || yield_stress < 1e-30 {
        return 0.0;
    }
    4.0 * k1 * k1 / (std::f64::consts::PI * young * yield_stress)
}
/// J-integral (energy release rate) for a mode-I crack in plane stress.
///
/// `J = K_I² / E`
pub fn j_integral_mode1(k1: f64, young: f64) -> f64 {
    if young < 1e-30 {
        return 0.0;
    }
    k1 * k1 / young
}
/// Crack branching angle by the maximum circumferential stress criterion.
///
/// Returns `θ_c` in radians for the direction of maximum σ_θθ at the crack tip.
///
/// From Williams (1957): `θ_c = 2 atan((K_I - sqrt(K_I^2 + 8 K_II^2)) / (4 K_II))`
/// when K_II ≠ 0, else 0.
pub fn branching_angle_max_circumferential(k1: f64, k2: f64) -> f64 {
    if k2.abs() < 1e-30 {
        return 0.0;
    }
    let disc = (k1 * k1 + 8.0 * k2 * k2).sqrt();
    2.0 * ((k1 - disc) / (4.0 * k2)).atan()
}
#[cfg(test)]
mod fracture_extended_tests {

    use crate::fracture::CohesiveZone;

    use crate::fracture::Fragment;

    use crate::fracture::branching_angle_max_circumferential;

    use crate::fracture::crack_opening_displacement;
    use crate::fracture::dynamic_fragment_count;

    use crate::fracture::extract_fragments;
    use crate::fracture::griffith_critical_stress;

    use crate::fracture::j_integral_mode1;

    use crate::fracture::mott_fragment_size;

    use crate::fracture::sif_edge_crack_finite_plate;

    #[test]
    fn test_cohesive_zone_zero_traction_at_zero_opening() {
        let mut cz = CohesiveZone::new(1e6, 0.001, 0.01);
        let t = cz.traction(0.0);
        assert_eq!(t, 0.0, "zero opening should give zero traction");
    }
    #[test]
    fn test_cohesive_zone_peak_at_delta0() {
        let mut cz = CohesiveZone::new(1e6, 0.001, 0.01);
        let t = cz.traction(0.001);
        assert!(
            (t - 1e6).abs() < 1e6 * 1e-10,
            "traction at δ₀ should equal T_c, got {t}"
        );
    }
    #[test]
    fn test_cohesive_zone_softening() {
        let mut cz = CohesiveZone::new(1e6, 0.001, 0.01);
        let t_peak = cz.traction(0.001);
        let _t_mid = cz.traction(0.005);
        let mut cz2 = CohesiveZone::new(1e6, 0.001, 0.01);
        cz2.traction(0.001);
        let t_soft = cz2.traction(0.005);
        assert!(
            t_soft < t_peak,
            "traction should decrease in softening branch"
        );
    }
    #[test]
    fn test_cohesive_zone_zero_after_separation() {
        let mut cz = CohesiveZone::new(1e6, 0.001, 0.01);
        cz.traction(0.015);
        assert!(cz.is_separated(), "should be separated after δ > δ_c");
        let t = cz.traction(0.001);
        assert_eq!(t, 0.0, "traction should be zero after full separation");
    }
    #[test]
    fn test_cohesive_zone_fracture_energy_formula() {
        let t_c = 1e6_f64;
        let delta_c = 0.01_f64;
        let cz = CohesiveZone::new(t_c, 0.001, delta_c);
        let g_c = cz.fracture_energy();
        let expected = 0.5 * t_c * delta_c;
        assert!(
            (g_c - expected).abs() < 1e-6,
            "G_c = 0.5 T_c δ_c, expected {expected}, got {g_c}"
        );
    }
    #[test]
    fn test_cohesive_zone_reset() {
        let mut cz = CohesiveZone::new(1e6, 0.001, 0.01);
        cz.traction(0.005);
        assert!(cz.delta_max > 0.0);
        cz.reset();
        assert_eq!(cz.delta_max, 0.0);
        assert!(!cz.separated);
    }
    #[test]
    fn test_cohesive_zone_compression_no_traction() {
        let mut cz = CohesiveZone::new(1e6, 0.001, 0.01);
        let t = cz.traction(-0.001);
        assert_eq!(t, 0.0, "compression should give zero traction");
    }
    #[test]
    fn test_cohesive_zone_secant_stiffness_linear_regime() {
        let cz = CohesiveZone::new(1e6, 0.001, 0.01);
        let k = cz.secant_stiffness(0.0005);
        let expected = 1e6 / 0.001;
        assert!(
            (k - expected).abs() / expected < 1e-6,
            "secant stiffness in linear regime, expected {expected}, got {k}"
        );
    }
    #[test]
    fn test_dynamic_fragment_count_minimum_one() {
        let n = dynamic_fragment_count(0.0, 1e6, 200e9, 0.01);
        assert_eq!(n, 1, "zero impact should give 1 fragment");
    }
    #[test]
    fn test_dynamic_fragment_count_increases_with_energy() {
        let n1 = dynamic_fragment_count(10.0, 1e6, 200e9, 0.01);
        let n2 = dynamic_fragment_count(100.0, 1e6, 200e9, 0.01);
        assert!(
            n2 >= n1,
            "more energy should give at least as many fragments"
        );
    }
    #[test]
    fn test_mott_fragment_size_decreases_with_strain_rate() {
        let l1 = mott_fragment_size(1e6, 7800.0, 200e9, 100.0);
        let l2 = mott_fragment_size(1e6, 7800.0, 200e9, 1000.0);
        assert!(
            l2 < l1,
            "higher strain rate should give smaller fragments, l1={l1}, l2={l2}"
        );
    }
    #[test]
    fn test_mott_fragment_size_increases_with_toughness() {
        let l1 = mott_fragment_size(1e6, 7800.0, 200e9, 1000.0);
        let l2 = mott_fragment_size(2e6, 7800.0, 200e9, 1000.0);
        assert!(l2 > l1, "higher toughness should give larger fragments");
    }
    #[test]
    fn test_extract_fragments_fully_connected() {
        let positions = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 3];
        let masses = [1.0; 3];
        let edges = [(0, 1), (1, 2)];
        let broken = [false, false];
        let frags = extract_fragments(&positions, &velocities, &masses, &edges, &broken);
        assert_eq!(frags.len(), 1, "fully connected → 1 fragment");
        assert_eq!(frags[0].node_indices.len(), 3);
    }
    #[test]
    fn test_extract_fragments_all_broken() {
        let positions = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 3];
        let masses = [1.0; 3];
        let edges = [(0, 1), (1, 2)];
        let broken = [true, true];
        let frags = extract_fragments(&positions, &velocities, &masses, &edges, &broken);
        assert_eq!(frags.len(), 3, "all broken → 3 fragments");
    }
    #[test]
    fn test_extract_fragments_partial_break() {
        let positions = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 4];
        let masses = [1.0; 4];
        let edges = [(0, 1), (1, 2), (2, 3)];
        let broken = [false, true, false];
        let frags = extract_fragments(&positions, &velocities, &masses, &edges, &broken);
        assert_eq!(frags.len(), 2, "one break → 2 fragments");
    }
    #[test]
    fn test_fragment_center_of_mass_midpoint() {
        let positions = [[0.0; 3], [2.0, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 2];
        let masses = [1.0, 1.0];
        let frag = Fragment::from_nodes(vec![0, 1], &positions, &velocities, &masses);
        assert!(
            (frag.center_of_mass[0] - 1.0).abs() < 1e-10,
            "COM should be at midpoint"
        );
    }
    #[test]
    fn test_fragment_kinetic_energy() {
        let positions = [[0.0; 3]];
        let velocities = [[3.0, 4.0, 0.0]];
        let masses = [2.0];
        let frag = Fragment::from_nodes(vec![0], &positions, &velocities, &masses);
        let ke = frag.kinetic_energy();
        assert!((ke - 25.0).abs() < 1e-10, "KE = 0.5*2*25 = 25, got {ke}");
    }
    #[test]
    fn test_sif_edge_crack_reduces_to_infinite_at_small_ratio() {
        let sigma = 100e6_f64;
        let a = 0.001_f64;
        let w = 1.0_f64;
        let k = sif_edge_crack_finite_plate(sigma, a, w);
        let k_infty = sigma * (std::f64::consts::PI * a).sqrt();
        assert!(
            (k / k_infty - 1.12).abs() < 0.05,
            "edge crack SIF should be ~1.12 × infinite-plate SIF for small a/W, got ratio {}",
            k / k_infty
        );
    }
    #[test]
    fn test_griffith_critical_stress_positive() {
        let sigma_c = griffith_critical_stress(200e9, 10.0, 0.01);
        assert!(
            sigma_c > 0.0,
            "critical stress should be positive, got {sigma_c}"
        );
    }
    #[test]
    fn test_griffith_critical_stress_decreases_with_crack_length() {
        let s1 = griffith_critical_stress(200e9, 10.0, 0.001);
        let s2 = griffith_critical_stress(200e9, 10.0, 0.01);
        assert!(s2 < s1, "critical stress decreases with crack length");
    }
    #[test]
    fn test_crack_opening_displacement_dugdale() {
        let k1 = 10e6_f64;
        let e = 200e9_f64;
        let sy = 300e6_f64;
        let cod = crack_opening_displacement(k1, e, sy);
        assert!(cod > 0.0, "COD should be positive");
        let expected = 4.0 * k1 * k1 / (std::f64::consts::PI * e * sy);
        assert!(
            (cod - expected).abs() < 1e-20,
            "COD formula check: {cod} vs {expected}"
        );
    }
    #[test]
    fn test_j_integral_mode1() {
        let k1 = 10e6_f64;
        let e = 200e9_f64;
        let j = j_integral_mode1(k1, e);
        let expected = k1 * k1 / e;
        assert!(
            (j - expected).abs() < 1e-15,
            "J = K²/E, expected {expected}, got {j}"
        );
    }
    #[test]
    fn test_branching_angle_pure_mode1_is_zero() {
        let theta = branching_angle_max_circumferential(1e6, 0.0);
        assert!(
            theta.abs() < 1e-12,
            "pure mode-I branching angle should be 0, got {theta}"
        );
    }
    #[test]
    fn test_branching_angle_nonzero_for_mixed_mode() {
        let theta = branching_angle_max_circumferential(1e6, 3e5);
        assert!(
            theta.abs() > 0.01,
            "mixed mode should give nonzero branching angle, got {theta}"
        );
    }
    #[test]
    fn test_branching_angle_antisymmetric_in_k2() {
        let theta_pos = branching_angle_max_circumferential(1e6, 3e5);
        let theta_neg = branching_angle_max_circumferential(1e6, -3e5);
        assert!(
            (theta_pos + theta_neg).abs() < 1e-10,
            "angle should be antisymmetric in K_II"
        );
    }
}
