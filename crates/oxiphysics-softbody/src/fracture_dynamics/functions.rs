//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CrackTip, FractureMaterial, FractureMesh};

pub(super) fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn vec_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn vec_length(a: [f64; 3]) -> f64 {
    vec_dot(a, a).sqrt()
}
pub(super) fn vec_normalize(a: [f64; 3]) -> [f64; 3] {
    let len = vec_length(a);
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        vec_scale(a, 1.0 / len)
    }
}
/// Build a rectangular tension specimen with a pre-existing edge notch on the
/// left side.
///
/// The specimen has width `width` (X direction) and height `height` (Y
/// direction).  A notch of depth `notch_depth` is cut from the left edge at
/// mid-height by removing / leaving-broken all bonds whose midpoints lie
/// within the notch region.  `resolution` controls the node spacing.
pub fn notched_tension_specimen(
    width: f64,
    height: f64,
    notch_depth: f64,
    resolution: f64,
    mat: FractureMaterial,
) -> FractureMesh {
    let dx = resolution;
    let nx = ((width / dx).round() as usize).max(2);
    let ny = ((height / dx).round() as usize).max(2);
    let mut mesh = FractureMesh::new_grid(nx, ny, dx, mat);
    let notch_y_lo = height / 2.0 - dx * 0.5;
    let notch_y_hi = height / 2.0 + dx * 0.5;
    let positions: Vec<[f64; 3]> = mesh.nodes.iter().map(|n| n.pos).collect();
    for bond in &mut mesh.bonds {
        let pa = positions[bond.node_a];
        let pb = positions[bond.node_b];
        let mid_x = (pa[0] + pb[0]) * 0.5;
        let mid_y = (pa[1] + pb[1]) * 0.5;
        if mid_x <= notch_depth && mid_y >= notch_y_lo && mid_y <= notch_y_hi {
            bond.broken = true;
        }
    }
    let notch_front_x = mesh
        .bonds
        .iter()
        .filter(|b| b.broken)
        .map(|b| {
            let pa = mesh.nodes[b.node_a].pos;
            let pb = mesh.nodes[b.node_b].pos;
            (pa[0] + pb[0]) * 0.5
        })
        .fold(f64::NEG_INFINITY, f64::max);
    if notch_front_x.is_finite() {
        let mid_y = height / 2.0;
        let tip_node = mesh
            .nodes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = (a.pos[0] - notch_front_x).powi(2) + (a.pos[1] - mid_y).powi(2);
                let db = (b.pos[0] - notch_front_x).powi(2) + (b.pos[1] - mid_y).powi(2);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        mesh.crack_tips
            .push(CrackTip::new(tip_node, [1.0, 0.0, 0.0]));
    }
    mesh
}
/// Compute the Rayleigh wave speed for a linear-elastic material.
///
/// Uses the Viktorov approximation: `c_R ≈ c_S * (0.862 + 1.14 ν) / (1 + ν)`
/// where `c_S = sqrt(G / ρ)` is the shear-wave speed.
pub fn rayleigh_wave_speed(mat: &FractureMaterial) -> f64 {
    let nu = mat.poisson_ratio;
    let e = mat.youngs_modulus;
    let rho = mat.density;
    let g_shear = e / (2.0 * (1.0 + nu));
    let c_s = (g_shear / rho).sqrt();
    c_s * (0.862 + 1.14 * nu) / (1.0 + nu)
}
/// Compute the longitudinal (P-wave) speed: `c_L = sqrt(E/ρ)` (plane stress).
pub fn longitudinal_wave_speed(mat: &FractureMaterial) -> f64 {
    (mat.youngs_modulus / mat.density).sqrt()
}
/// Compute the dynamic stress intensity factor (SIF) for a crack tip moving
/// at velocity `v` relative to the Rayleigh wave speed `c_R`.
///
/// Uses the Freund formula:
/// `K_I^dyn = K_I^static * k(v) ≈ K_I^static * (1 - v/c_R)`
pub fn dynamic_stress_intensity(k_i_static: f64, crack_speed: f64, c_rayleigh: f64) -> f64 {
    let ratio = (crack_speed / c_rayleigh).min(0.99);
    k_i_static * (1.0 - ratio)
}
/// Estimate the crack propagation speed from the applied SIF and fracture
/// toughness using a simplified power law: `v = c_R * (1 - K_Ic / K_I)`.
///
/// Returns 0 when `K_I <= K_Ic`.
pub fn crack_velocity(k_i: f64, k_ic: f64, c_rayleigh: f64) -> f64 {
    if k_i <= k_ic {
        0.0
    } else {
        c_rayleigh * (1.0 - k_ic / k_i).max(0.0)
    }
}
/// Branching criterion: crack branches when the crack velocity exceeds the
/// branching threshold, typically `≈ 0.4 * c_R` (Yoffe criterion).
pub fn should_branch(crack_speed: f64, c_rayleigh: f64, threshold_fraction: f64) -> bool {
    crack_speed >= threshold_fraction * c_rayleigh
}
/// Generate two branching directions at ±`branch_angle_deg` degrees from the
/// current crack direction.
#[allow(clippy::too_many_arguments)]
pub fn branch_directions(current_dir: [f64; 3], branch_angle_deg: f64) -> ([f64; 3], [f64; 3]) {
    let theta = branch_angle_deg.to_radians();
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let dx = current_dir[0];
    let dy = current_dir[1];
    let dir_pos = vec_normalize([
        dx * cos_t - dy * sin_t,
        dx * sin_t + dy * cos_t,
        current_dir[2],
    ]);
    let dir_neg = vec_normalize([
        dx * cos_t + dy * sin_t,
        -dx * sin_t + dy * cos_t,
        current_dir[2],
    ]);
    (dir_pos, dir_neg)
}
#[cfg(test)]
mod tests {
    use super::*;

    fn glass_material() -> FractureMaterial {
        FractureMaterial::new(70e9, 0.23, 2500.0, 50e6, 0.75)
    }
    #[test]
    fn test_new_grid_node_count() {
        let mat = glass_material();
        let mesh = FractureMesh::new_grid(5, 4, 0.1, mat);
        assert_eq!(mesh.nodes.len(), 20, "5×4 grid should have 20 nodes");
    }
    #[test]
    fn test_compute_bond_forces_undeformed() {
        let mat = glass_material();
        let mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        let forces = mesh.compute_bond_forces();
        for (i, f) in forces.iter().enumerate() {
            let mag = vec_length(*f);
            assert!(
                mag < 1e-10,
                "Node {i}: force {mag} should be zero for undeformed mesh"
            );
        }
    }
    #[test]
    fn test_check_fracture_breaks_overstretched_bond() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(2, 2, 0.1, mat);
        mesh.nodes[1].pos[0] = 10.0;
        mesh.check_fracture();
        let broken = mesh.bonds.iter().filter(|b| b.broken).count();
        assert!(
            broken > 0,
            "At least one bond should break after extreme stretch"
        );
    }
    #[test]
    fn test_broken_bond_count_increases() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(2, 2, 0.1, mat);
        let before = mesh.broken_bond_count();
        mesh.nodes[3].pos = [100.0, 100.0, 0.0];
        mesh.check_fracture();
        let after = mesh.broken_bond_count();
        assert!(
            after > before,
            "broken_bond_count should increase after fracture: before={before}, after={after}"
        );
    }
    #[test]
    fn test_integrate_nodes_move_under_gravity() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(2, 2, 0.1, mat);
        let initial_y: Vec<f64> = mesh.nodes.iter().map(|n| n.pos[1]).collect();
        let dt = 1.0 / 60.0;
        for _ in 0..10 {
            mesh.integrate(dt);
        }
        let final_y: Vec<f64> = mesh.nodes.iter().map(|n| n.pos[1]).collect();
        let moved = initial_y
            .iter()
            .zip(final_y.iter())
            .any(|(y0, y1)| (y1 - y0).abs() > 1e-10);
        assert!(moved, "Nodes should move downward under gravity");
    }
    #[test]
    fn test_damage_field_after_fracture() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        mesh.nodes[1].pos[0] = 50.0;
        mesh.check_fracture();
        let damages = mesh.damage_field();
        let any_damaged = damages.iter().any(|&d| d > 0.0);
        assert!(
            any_damaged,
            "At least one node should have non-zero damage after fracture"
        );
    }
    #[test]
    fn test_total_fracture_energy_non_negative() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        mesh.nodes[4].pos = [1000.0, 1000.0, 0.0];
        mesh.check_fracture();
        let e = mesh.total_fracture_energy();
        assert!(
            e >= 0.0,
            "Total fracture energy must be non-negative, got {e}"
        );
    }
    #[test]
    fn test_rayleigh_wave_speed() {
        let mat = glass_material();
        let c_r = rayleigh_wave_speed(&mat);
        let c_l = longitudinal_wave_speed(&mat);
        assert!(c_r > 0.0, "Rayleigh speed must be positive, got {c_r}");
        assert!(
            c_r < c_l,
            "Rayleigh speed {c_r} must be < longitudinal speed {c_l}"
        );
    }
    #[test]
    fn test_dynamic_sif_decreases_with_speed() {
        let mat = glass_material();
        let c_r = rayleigh_wave_speed(&mat);
        let k_static = 1.0e6;
        let k_slow = dynamic_stress_intensity(k_static, 0.1 * c_r, c_r);
        let k_fast = dynamic_stress_intensity(k_static, 0.5 * c_r, c_r);
        assert!(
            k_slow > k_fast,
            "Dynamic SIF should decrease with speed: {k_slow} vs {k_fast}"
        );
    }
    #[test]
    fn test_crack_velocity_below_kic() {
        let mat = glass_material();
        let c_r = rayleigh_wave_speed(&mat);
        let v = crack_velocity(0.5 * mat.fracture_toughness, mat.fracture_toughness, c_r);
        assert_eq!(v, 0.0, "crack velocity must be zero when K_I < K_Ic");
    }
    #[test]
    fn test_crack_velocity_above_kic() {
        let mat = glass_material();
        let c_r = rayleigh_wave_speed(&mat);
        let v = crack_velocity(2.0 * mat.fracture_toughness, mat.fracture_toughness, c_r);
        assert!(
            v > 0.0,
            "crack velocity must be positive when K_I > K_Ic, got {v}"
        );
    }
    #[test]
    fn test_branch_directions_unit_length() {
        let dir = [1.0, 0.0, 0.0];
        let (d1, d2) = branch_directions(dir, 30.0);
        let l1 = vec_length(d1);
        let l2 = vec_length(d2);
        assert!((l1 - 1.0).abs() < 1e-12, "Branch dir 1 not unit: {l1}");
        assert!((l2 - 1.0).abs() < 1e-12, "Branch dir 2 not unit: {l2}");
    }
    #[test]
    fn test_branch_directions_angle() {
        let dir = [1.0, 0.0, 0.0];
        let angle = 30.0_f64;
        let (d1, _d2) = branch_directions(dir, angle);
        let cos_angle = vec_dot(dir, d1);
        let expected_cos = angle.to_radians().cos();
        assert!(
            (cos_angle - expected_cos).abs() < 1e-10,
            "Branch angle incorrect: cos={cos_angle}, expected {expected_cos}"
        );
    }
    #[test]
    fn test_should_branch() {
        let mat = glass_material();
        let c_r = rayleigh_wave_speed(&mat);
        assert!(!should_branch(0.3 * c_r, c_r, 0.4));
        assert!(should_branch(0.5 * c_r, c_r, 0.4));
    }
    #[test]
    fn test_try_branch_cracks_spawns_tips() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(4, 4, 0.1, mat);
        let c_r = rayleigh_wave_speed(&mesh.material);
        let mut tip = CrackTip::new(5, [1.0, 0.0, 0.0]);
        tip.velocity = vec_scale([1.0, 0.0, 0.0], 0.5 * c_r);
        mesh.crack_tips.push(tip);
        let before = mesh.crack_tips.len();
        mesh.try_branch_cracks(30.0, 0.4);
        let after = mesh.crack_tips.len();
        assert!(
            after > before,
            "Branch should spawn new crack tips: before={before}, after={after}"
        );
    }
    #[test]
    fn test_total_strain_energy() {
        let mat = glass_material();
        let mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        let e = mesh.total_strain_energy();
        assert!(e >= 0.0, "Strain energy must be non-negative, got {e}");
    }
    #[test]
    fn test_kinetic_energy() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        for _ in 0..5 {
            mesh.integrate(1e-4);
        }
        let ke = mesh.kinetic_energy();
        assert!(ke >= 0.0, "Kinetic energy must be non-negative, got {ke}");
    }
    #[test]
    fn test_apply_damping_reduces_ke() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        for _ in 0..10 {
            mesh.integrate(1e-4);
        }
        let ke_before = mesh.kinetic_energy();
        mesh.apply_damping(0.5);
        let ke_after = mesh.kinetic_energy();
        assert!(
            ke_after <= ke_before + 1e-20,
            "Damping should reduce KE: {ke_before} → {ke_after}"
        );
    }
    #[test]
    fn test_center_of_mass_in_bounds() {
        let mat = glass_material();
        let dx = 0.1;
        let mesh = FractureMesh::new_grid(4, 4, dx, mat);
        let com = mesh.center_of_mass();
        let xlo = mesh
            .nodes
            .iter()
            .map(|n| n.pos[0])
            .fold(f64::INFINITY, f64::min);
        let xhi = mesh
            .nodes
            .iter()
            .map(|n| n.pos[0])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            com[0] >= xlo - 1e-10 && com[0] <= xhi + 1e-10,
            "CoM x={} out of bounds [{xlo}, {xhi}]",
            com[0]
        );
    }
    #[test]
    fn test_apply_tensile_load_y() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        mesh.apply_tensile_load_y(0.01);
        let fixed_count = mesh.nodes.iter().filter(|n| n.fixed).count();
        assert!(fixed_count > 0, "Tensile load should fix some nodes");
    }
    #[test]
    fn test_notched_specimen_has_crack_tip() {
        let mat = glass_material();
        let specimen = notched_tension_specimen(1.0, 0.5, 0.2, 0.1, mat);
        assert!(
            !specimen.crack_tips.is_empty(),
            "Notched specimen should have at least one crack tip"
        );
    }
    #[test]
    fn test_bond_count_consistency() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        mesh.nodes[4].pos = [10.0, 10.0, 0.0];
        mesh.check_fracture();
        assert_eq!(
            mesh.intact_bond_count() + mesh.broken_bond_count(),
            mesh.bonds.len()
        );
    }
    #[test]
    fn test_average_damage_range() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(4, 4, 0.1, mat);
        mesh.nodes[5].pos = [500.0, 500.0, 0.0];
        mesh.check_fracture();
        let d = mesh.average_damage();
        assert!(
            (0.0..=1.0).contains(&d),
            "average_damage must be in [0,1], got {d}"
        );
    }
    #[test]
    fn test_step_substeps_no_panic() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(3, 3, 0.1, mat);
        mesh.step_substeps(1e-3, 4);
        assert!(mesh.kinetic_energy().is_finite());
    }
    #[test]
    fn test_apply_nodal_force() {
        let mat = glass_material();
        let mut mesh = FractureMesh::new_grid(2, 2, 0.1, mat);
        let vel_before = mesh.nodes[0].vel;
        let force = [1.0e6, 0.0, 0.0];
        mesh.apply_nodal_force(0, force, 1e-4);
        let vel_after = mesh.nodes[0].vel;
        assert!(
            vel_after[0] > vel_before[0],
            "Force should increase x-velocity"
        );
    }
}
/// Evaluate the Griffith criterion for a crack of length `a`.
///
/// The Griffith criterion states that a crack of half-length `a` will propagate
/// when the stress intensity factor K_I ≥ K_Ic:
///
/// ```text
/// K_I = sigma * sqrt(pi * a)
/// ```
///
/// Returns `true` if the crack will propagate.
#[allow(dead_code)]
pub fn griffith_criterion(sigma: f64, crack_half_length: f64, k_ic: f64) -> bool {
    if crack_half_length <= 0.0 || sigma <= 0.0 {
        return false;
    }
    let k_i = sigma * (std::f64::consts::PI * crack_half_length).sqrt();
    k_i >= k_ic
}
/// Compute the critical half-crack length at which fracture initiates under
/// applied stress `sigma` (Griffith critical crack size).
///
/// `a_c = (K_Ic / (sigma * sqrt(pi)))^2`
#[allow(dead_code)]
pub fn griffith_critical_crack_size(sigma: f64, k_ic: f64) -> f64 {
    if sigma < 1e-15 {
        return f64::INFINITY;
    }
    let ratio = k_ic / (sigma * std::f64::consts::PI.sqrt());
    ratio * ratio
}
/// Compute the energy release rate G for a plane-stress crack.
///
/// `G = K_I^2 / E` (plane stress)
/// `G = K_I^2 * (1 - nu^2) / E` (plane strain)
#[allow(dead_code)]
pub fn energy_release_rate(
    k_i: f64,
    youngs_modulus: f64,
    poisson_ratio: f64,
    plane_strain: bool,
) -> f64 {
    if youngs_modulus < 1e-15 {
        return 0.0;
    }
    if plane_strain {
        k_i * k_i * (1.0 - poisson_ratio * poisson_ratio) / youngs_modulus
    } else {
        k_i * k_i / youngs_modulus
    }
}
/// Compute crack branching probability based on the dynamic stress intensity factor.
///
/// Branching occurs when the crack velocity exceeds a fraction of the Rayleigh wave
/// speed cR. The branching probability increases from 0 at v = 0.4*cR to 1 at v = cR.
///
/// `P_branch = clamp((v_crack/cR - 0.4) / 0.6, 0, 1)`
#[allow(dead_code)]
pub fn branching_probability(crack_velocity: f64, rayleigh_wave_speed: f64) -> f64 {
    if rayleigh_wave_speed < 1e-15 {
        return 0.0;
    }
    let v_norm = crack_velocity / rayleigh_wave_speed;
    ((v_norm - 0.4) / 0.6).clamp(0.0, 1.0)
}
/// Compute the Rayleigh wave speed from elastic constants.
///
/// The Rayleigh wave speed cR satisfies:
/// `cR ≈ cs * (0.862 + 1.14*nu) / (1 + nu)`
///
/// where cs = sqrt(G/rho) is the shear wave speed and G = E / (2*(1+nu)).
#[allow(dead_code)]
pub fn rayleigh_wave_speed_scalar(youngs_modulus: f64, poisson_ratio: f64, density: f64) -> f64 {
    if density < 1e-15 || youngs_modulus < 1e-15 {
        return 0.0;
    }
    let nu = poisson_ratio;
    let g = youngs_modulus / (2.0 * (1.0 + nu));
    let cs = (g / density).sqrt();
    cs * (0.862 + 1.14 * nu) / (1.0 + nu)
}
/// Stochastically determine whether branching occurs based on probability and a threshold.
///
/// Returns `true` if a random sample `r` ∈ \[0,1) is less than `prob`.
#[allow(dead_code)]
pub fn sample_branching(prob: f64, random_sample: f64) -> bool {
    random_sample < prob
}
#[cfg(test)]
mod tests_energy {

    use crate::CodProfile;
    use crate::CohesiveZone;
    use crate::CrackBranching;
    use crate::CrackPropagation;
    use crate::DynamicEnergyBalance;
    use crate::DynamicFractureToughness;
    use crate::FractureBond;
    use crate::FractureNode;
    use crate::FragmentDistribution;
    use crate::branching_probability;
    use crate::energy_release_rate;
    use crate::griffith_criterion;
    use crate::griffith_critical_crack_size;
    use crate::rayleigh_wave_speed_scalar;
    use crate::sample_branching;
    #[test]
    fn test_energy_balance_kinetic() {
        let mut eb = DynamicEnergyBalance::new();
        let nodes = vec![
            FractureNode {
                pos: [0.0; 3],
                vel: [1.0, 0.0, 0.0],
                mass: 2.0,
                fixed: false,
                damage: 0.0,
            },
            FractureNode {
                pos: [1.0, 0.0, 0.0],
                vel: [0.0, 2.0, 0.0],
                mass: 1.0,
                fixed: false,
                damage: 0.0,
            },
        ];
        eb.update_kinetic(&nodes);
        assert!(
            (eb.kinetic_energy - 3.0).abs() < 1e-12,
            "KE = {}",
            eb.kinetic_energy
        );
    }
    #[test]
    fn test_energy_balance_strain() {
        let mut eb = DynamicEnergyBalance::new();
        let bonds = vec![FractureBond {
            node_a: 0,
            node_b: 1,
            rest_length: 1.0,
            stiffness: 1000.0,
            strength: 1e6,
            fracture_energy: 100.0,
            broken: false,
            current_stretch: 0.1,
        }];
        eb.update_strain(&bonds);
        assert!(
            (eb.strain_energy - 5.0).abs() < 1e-10,
            "SE = {}",
            eb.strain_energy
        );
    }
    #[test]
    fn test_energy_balance_broken_bond_excluded() {
        let mut eb = DynamicEnergyBalance::new();
        let bonds = vec![FractureBond {
            node_a: 0,
            node_b: 1,
            rest_length: 1.0,
            stiffness: 1000.0,
            strength: 1e6,
            fracture_energy: 100.0,
            broken: true,
            current_stretch: 0.5,
        }];
        eb.update_strain(&bonds);
        assert!(
            (eb.strain_energy).abs() < 1e-14,
            "Broken bond should not contribute"
        );
    }
    #[test]
    fn test_energy_balance_record_and_drift() {
        let mut eb = DynamicEnergyBalance::new();
        eb.kinetic_energy = 50.0;
        eb.strain_energy = 50.0;
        eb.record();
        eb.kinetic_energy = 48.0;
        eb.strain_energy = 52.0;
        eb.record();
        let drift = eb.relative_drift();
        assert!(drift < 1e-12, "Conservative system: drift = {drift}");
    }
    #[test]
    fn test_energy_balance_fracture_increases_total() {
        let mut eb = DynamicEnergyBalance::new();
        eb.kinetic_energy = 100.0;
        eb.strain_energy = 50.0;
        eb.record();
        eb.add_fracture_event(10.0, 0.01);
        eb.kinetic_energy = 88.0;
        eb.strain_energy = 49.0;
        eb.record();
        assert!(eb.total_energy() > 0.0);
    }
    #[test]
    fn test_griffith_criterion() {
        let k_ic = 0.7e6;
        let sigma = 50e6;
        let a_c = griffith_critical_crack_size(sigma, k_ic);
        assert!(
            griffith_criterion(sigma, a_c, k_ic),
            "at critical size should propagate"
        );
        assert!(
            !griffith_criterion(sigma, a_c * 0.5, k_ic),
            "below critical should not propagate"
        );
    }
    #[test]
    fn test_griffith_zero_stress() {
        assert!(!griffith_criterion(0.0, 0.01, 1e6));
    }
    #[test]
    fn test_energy_release_rate_plane_stress() {
        let e = 70e9;
        let k_i = 1e6;
        let g = energy_release_rate(k_i, e, 0.33, false);
        let expected = k_i * k_i / e;
        assert!((g - expected).abs() < 1e-20, "G = {g}, expected {expected}");
    }
    #[test]
    fn test_energy_release_rate_plane_strain() {
        let e = 70e9;
        let nu = 0.33;
        let k_i = 1e6;
        let g_ps = energy_release_rate(k_i, e, nu, false);
        let g_pe = energy_release_rate(k_i, e, nu, true);
        assert!(
            g_pe < g_ps,
            "Plane strain G should be less than plane stress G"
        );
    }
    #[test]
    fn test_cohesive_zone_traction_linear() {
        let cz = CohesiveZone::new(100.0e6, 1e-4);
        let t0 = cz.traction(0.0);
        let t_half = cz.traction(0.5e-4);
        let t_full = cz.traction(1e-4);
        assert!(t0.abs() < 1e-10, "T(0) should be 0");
        assert!(
            (t_half - 50.0e6).abs() < 1.0,
            "T(delta_c/2) should be T_max/2"
        );
        assert!(t_full.abs() < 1e-10, "T(delta_c) should be 0");
    }
    #[test]
    fn test_cohesive_zone_update_failure() {
        let mut cz = CohesiveZone::new(100.0e6, 1e-4);
        let t = cz.update(2e-4);
        assert!(t.abs() < 1e-10, "Should be failed");
        assert!(cz.failed);
    }
    #[test]
    fn test_cohesive_zone_damage() {
        let mut cz = CohesiveZone::new(100.0e6, 1e-4);
        cz.update(5e-5);
        let d = cz.damage();
        assert!((d - 0.5).abs() < 1e-10, "Damage should be 0.5, got {d}");
    }
    #[test]
    fn test_cohesive_zone_fracture_energy() {
        let cz = CohesiveZone::new(100.0, 1e-3);
        let g = cz.fracture_energy();
        let expected = 0.5 * 100.0 * 1e-3;
        assert!((g - expected).abs() < 1e-15, "Gc = {g}");
    }
    #[test]
    fn test_cod_profile_interpolation() {
        let mut cod = CodProfile::new();
        cod.add_point(0.0, 0.0);
        cod.add_point(0.5, 0.5e-3);
        cod.add_point(1.0, 1e-3);
        let v = cod.interpolate(0.25);
        assert!((v - 0.25e-3).abs() < 1e-15, "Interpolated = {v}");
    }
    #[test]
    fn test_cod_profile_ctod() {
        let mut cod = CodProfile::new();
        cod.add_point(0.0, 1e-6);
        cod.add_point(0.5, 5e-4);
        cod.add_point(1.0, 1e-3);
        assert!((cod.ctod() - 1e-6).abs() < 1e-20, "CTOD should be at tip");
        assert!((cod.max_opening() - 1e-3).abs() < 1e-20, "Max at mouth");
    }
    #[test]
    fn test_branching_probability() {
        let cr = 3000.0;
        assert!(
            (branching_probability(0.0, cr)).abs() < 1e-10,
            "No motion = no branching"
        );
        assert!(
            (branching_probability(0.4 * cr, cr)).abs() < 1e-10,
            "Below threshold = 0"
        );
        let p_high = branching_probability(cr, cr);
        assert!((p_high - 1.0).abs() < 1e-6, "At cR = 1.0, got {p_high}");
    }
    #[test]
    fn test_rayleigh_wave_speed() {
        let e = 70e9;
        let nu = 0.33;
        let rho = 2700.0;
        let cr = rayleigh_wave_speed_scalar(e, nu, rho);
        assert!(cr > 0.0 && cr < 1e6, "Reasonable cR: {cr}");
    }
    #[test]
    fn test_dynamic_fracture_toughness_static_limit() {
        let kd = DynamicFractureToughness::new(1.0e6, 3000.0, 0.5);
        let k_id_zero = kd.k_id(0.0);
        assert!((k_id_zero - 1.0e6).abs() < 1e-6, "K_Id(0) = K_Ic");
    }
    #[test]
    fn test_dynamic_fracture_toughness_decreases() {
        let kd = DynamicFractureToughness::new(1.0e6, 3000.0, 0.5);
        let k_slow = kd.k_id(100.0);
        let k_fast = kd.k_id(2000.0);
        assert!(k_fast < k_slow, "K_Id should decrease with velocity");
    }
    #[test]
    fn test_sample_branching() {
        assert!(sample_branching(0.9, 0.5), "Should branch");
        assert!(!sample_branching(0.1, 0.5), "Should not branch");
        assert!(!sample_branching(0.0, 0.0), "Zero prob = no branch");
    }
    #[test]
    fn test_crack_propagation_below_kic_returns_zero() {
        let cp = CrackPropagation::new(1.0e6, 3000.0, 0.5, 0.5e6);
        let v = cp.compute_crack_tip_velocity(50);
        assert_eq!(v, 0.0, "Crack tip velocity should be 0 below K_Ic, got {v}");
    }
    #[test]
    fn test_crack_propagation_above_kic_positive() {
        let cp = CrackPropagation::new(1.0e6, 3000.0, 0.5, 3.0e6);
        let v = cp.compute_crack_tip_velocity(50);
        assert!(
            v > 0.0,
            "Crack tip velocity should be positive above K_Ic, got {v}"
        );
    }
    #[test]
    fn test_crack_propagation_below_rayleigh() {
        let c_r = 3000.0;
        let cp = CrackPropagation::new(1.0e6, c_r, 0.5, 5.0e6);
        let v = cp.compute_crack_tip_velocity(100);
        assert!(
            v < c_r,
            "Crack tip cannot exceed Rayleigh speed: v={v}, cR={c_r}"
        );
    }
    #[test]
    fn test_branching_angle_pure_mode_i() {
        let cb = CrackBranching::new(1.0e6, 0.0);
        let theta = cb.compute_branching_angle();
        assert!(
            theta.abs() < 1e-12,
            "Pure mode-I should give zero branching angle, got {theta}"
        );
    }
    #[test]
    fn test_branching_angle_mixed_mode_nonzero() {
        let cb = CrackBranching::new(1.0e6, 0.5e6);
        let theta = cb.compute_branching_angle();
        assert!(
            theta.abs() > 1e-6,
            "Mixed mode should give nonzero branching angle, got {theta}"
        );
    }
    #[test]
    fn test_branching_angle_clamped() {
        let cb = CrackBranching::new(1.0e6, 10.0e6);
        let theta = cb.compute_branching_angle();
        assert!(
            theta.abs() <= std::f64::consts::FRAC_PI_2 + 1e-12,
            "Branching angle must be within ±90 deg, got {}",
            theta.to_degrees()
        );
    }
    #[test]
    fn test_mott_grady_positive_size() {
        let fd = FragmentDistribution::new(1.0e6, 2500.0, 70e9);
        let lf = fd.compute_mott_grady(1000.0);
        assert!(
            lf > 0.0,
            "Mott-Grady fragment size should be positive, got {lf}"
        );
    }
    #[test]
    fn test_mott_grady_zero_strain_rate() {
        let fd = FragmentDistribution::new(1.0e6, 2500.0, 70e9);
        let lf = fd.compute_mott_grady(0.0);
        assert_eq!(lf, 0.0, "Mott-Grady should return 0 at zero strain rate");
    }
    #[test]
    fn test_mott_grady_decreases_with_strain_rate() {
        let fd = FragmentDistribution::new(1.0e6, 2500.0, 70e9);
        let lf_slow = fd.compute_mott_grady(100.0);
        let lf_fast = fd.compute_mott_grady(10_000.0);
        assert!(
            lf_fast < lf_slow,
            "Faster strain rate → smaller fragments: slow={lf_slow}, fast={lf_fast}"
        );
    }
    #[test]
    fn test_fragment_number_density_increases_with_strain_rate() {
        let fd = FragmentDistribution::new(1.0e6, 2500.0, 70e9);
        let n_slow = fd.fragment_number_density(100.0);
        let n_fast = fd.fragment_number_density(10_000.0);
        assert!(
            n_fast > n_slow,
            "Faster strain rate → more fragments per volume: slow={n_slow}, fast={n_fast}"
        );
    }
}
