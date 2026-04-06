//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CgBead, CgWaterBead, MartiniLj};

/// Dot product of two 3-vectors.
#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Subtraction of two 3-vectors.
#[inline]
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Addition of two 3-vectors.
#[inline]
pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Euclidean length of a 3-vector.
#[inline]
pub fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
/// Normalize a 3-vector to unit length. Returns zero vector if length is zero.
#[inline]
pub fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / len)
    }
}
/// Compute LJ interaction energy between a MARTINI W bead pair.
///
/// Uses sigma = 0.47 nm, epsilon = 5.0 kJ/mol (MARTINI W-W).
#[allow(dead_code)]
pub fn cg_water_energy(r: f64) -> f64 {
    MartiniLj::energy(r, 0.47, 5.0)
}
/// Build a simple cubic CG water box with `n` beads per side.
///
/// Returns bead positions arranged on a cubic grid with spacing `spacing` nm.
#[allow(dead_code)]
pub fn build_cg_water_box(n: usize, spacing: f64) -> Vec<CgWaterBead> {
    let mut beads = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let pos = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                beads.push(CgWaterBead::new(pos));
            }
        }
    }
    beads
}
/// Reconstruct approximate all-atom positions from CG bead positions.
///
/// This simple scheme places N atomic positions around a CG bead centroid
/// on a sphere of radius `r_aa` nm.
#[allow(dead_code)]
pub fn backmap_bead(bead_pos: [f64; 3], n_atoms: usize, r_aa: f64) -> Vec<[f64; 3]> {
    if n_atoms == 0 {
        return Vec::new();
    }
    if n_atoms == 1 {
        return vec![bead_pos];
    }
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..n_atoms)
        .map(|i| {
            let y = 1.0 - (i as f64 / (n_atoms - 1).max(1) as f64) * 2.0;
            let radius = (1.0 - y * y).max(0.0).sqrt();
            let theta = golden_angle * i as f64;
            [
                bead_pos[0] + r_aa * radius * theta.cos(),
                bead_pos[1] + r_aa * y,
                bead_pos[2] + r_aa * radius * theta.sin(),
            ]
        })
        .collect()
}
/// Estimate all-atom positions from multiple CG beads, distributing `atoms_per_bead`
/// atoms around each bead centroid.
#[allow(dead_code)]
pub fn backmap_molecule(
    bead_positions: &[[f64; 3]],
    atoms_per_bead: usize,
    r_aa: f64,
) -> Vec<[f64; 3]> {
    let mut all_positions = Vec::new();
    for &bead_pos in bead_positions {
        let positions = backmap_bead(bead_pos, atoms_per_bead, r_aa);
        all_positions.extend(positions);
    }
    all_positions
}
/// Compute the radius of gyration of a set of CG beads.
#[allow(dead_code)]
pub fn cg_radius_of_gyration(beads: &[CgBead]) -> f64 {
    if beads.is_empty() {
        return 0.0;
    }
    let total_mass: f64 = beads.iter().map(|b| b.mass).sum();
    if total_mass <= 0.0 {
        return 0.0;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for b in beads {
        cx += b.mass * b.pos[0];
        cy += b.mass * b.pos[1];
        cz += b.mass * b.pos[2];
    }
    cx /= total_mass;
    cy /= total_mass;
    cz /= total_mass;
    let sum_mr2: f64 = beads
        .iter()
        .map(|b| {
            let dx = b.pos[0] - cx;
            let dy = b.pos[1] - cy;
            let dz = b.pos[2] - cz;
            b.mass * (dx * dx + dy * dy + dz * dz)
        })
        .sum();
    (sum_mr2 / total_mass).sqrt()
}
/// Compute the end-to-end distance of a linear CG polymer chain.
#[allow(dead_code)]
pub fn end_to_end_distance(beads: &[CgBead]) -> f64 {
    if beads.len() < 2 {
        return 0.0;
    }
    let first = &beads[0];
    let last = &beads[beads.len() - 1];
    length(sub(last.pos, first.pos))
}
/// Compute the total kinetic energy of a set of CG beads.
#[allow(dead_code)]
pub fn cg_kinetic_energy(beads: &[CgBead]) -> f64 {
    beads
        .iter()
        .map(|b| {
            let v2 = dot(b.vel, b.vel);
            0.5 * b.mass * v2
        })
        .sum()
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use std::f64::consts::PI;
    fn make_bead(pos: [f64; 3]) -> CgBead {
        CgBead::new("B", 0, 72.0, 0.0, pos)
    }
    #[test]
    fn test_cgbond_energy_at_equilibrium() {
        let beads = vec![make_bead([0.0, 0.0, 0.0]), make_bead([0.47, 0.0, 0.0])];
        let bond = CgBond::new(0, 1, 0.47, 1000.0);
        let e = bond.energy(&beads);
        assert!(
            e.abs() < 1e-12,
            "Energy at equilibrium should be zero, got {}",
            e
        );
    }
    #[test]
    fn test_cgbond_forces_stretched() {
        let r0 = 0.47;
        let r = 0.60;
        let beads = vec![make_bead([0.0, 0.0, 0.0]), make_bead([r, 0.0, 0.0])];
        let bond = CgBond::new(0, 1, r0, 1000.0);
        let (fi, fj) = bond.forces(&beads);
        assert!(
            fi[0] > 0.0,
            "Force on i should be toward j, got fi[0]={}",
            fi[0]
        );
        assert!(
            fj[0] < 0.0,
            "Force on j should be toward i, got fj[0]={}",
            fj[0]
        );
        assert!((fi[0] + fj[0]).abs() < 1e-12);
    }
    #[test]
    fn test_enm_contact_map_close_beads() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let enm = ElasticNetworkModel::new(positions, 0.5, 10.0);
        let contacts = enm.contact_map();
        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0], (0, 1));
    }
    #[test]
    fn test_enm_energy_zero_at_reference() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let enm = ElasticNetworkModel::new(positions.clone(), 0.5, 10.0);
        let e = enm.potential_energy(&positions);
        assert!(
            e.abs() < 1e-12,
            "Energy at reference should be zero, got {}",
            e
        );
    }
    #[test]
    fn test_cgsystem_linear_chain_bond_count() {
        let n = 5;
        let sys = PolymerChain::linear_chain(n, 72.0, 0.47, 1000.0);
        assert_eq!(sys.beads.len(), n);
        assert_eq!(sys.bonds.len(), n - 1);
    }
    #[test]
    fn test_cgsystem_ring_polymer_bond_count() {
        let n = 6;
        let sys = PolymerChain::ring_polymer(n, 72.0, 0.47, 1000.0);
        assert_eq!(sys.beads.len(), n);
        assert_eq!(sys.bonds.len(), n);
    }
    #[test]
    fn test_cgsystem_step_verlet_beads_move() {
        let mut sys = CgSystem::new(2.0);
        sys.add_bead(CgBead::new("B0", 0, 72.0, 0.0, [0.0, 0.0, 0.0]));
        sys.add_bead(CgBead::new("B1", 0, 72.0, 0.0, [0.60, 0.0, 0.0]));
        sys.bonds.push(CgBond::new(0, 1, 0.47, 1000.0));
        sys.compute_forces();
        let pos0_before = sys.beads[0].pos;
        let pos1_before = sys.beads[1].pos;
        sys.step_verlet(0.001);
        let moved = sys.beads[0].pos[0] != pos0_before[0] || sys.beads[1].pos[0] != pos1_before[0];
        assert!(moved, "Beads should move after Verlet step");
    }
    #[test]
    fn test_martini_lj_energy_minimum() {
        let sigma = 0.47;
        let epsilon = 5.6;
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let e_min = MartiniLj::energy(r_min, sigma, epsilon);
        assert!(
            (e_min + epsilon).abs() < 1e-10,
            "LJ minimum energy should be -epsilon, got {}",
            e_min
        );
        let e_sigma = MartiniLj::energy(sigma, sigma, epsilon);
        assert!(
            e_sigma.abs() < 1e-10,
            "LJ energy at r=sigma should be 0, got {}",
            e_sigma
        );
    }
    #[test]
    fn test_martini_lj_energy_negative_between_sigma_and_rmin() {
        let sigma = 0.47;
        let epsilon = 5.6;
        let r = sigma * 1.1;
        let e = MartiniLj::energy(r, sigma, epsilon);
        assert!(
            e < 0.0,
            "LJ energy between sigma and r_min should be negative, got {}",
            e
        );
    }
    #[test]
    fn test_dihedral_angle_known_geometry() {
        let beads = vec![
            make_bead([0.0, 1.0, 0.0]),
            make_bead([0.0, 0.0, 0.0]),
            make_bead([1.0, 0.0, 0.0]),
            make_bead([1.0, 1.0, 0.0]),
        ];
        let dih = CgDihedral::new(0, 1, 2, 3, 0.0, 1.0, 1);
        let phi = dih.dihedral_angle(&beads);
        assert!(
            phi.abs() < 1e-10 || (phi - PI).abs() < 1e-10 || (phi + PI).abs() < 1e-10,
            "Expected known dihedral angle, got {}",
            phi
        );
    }
    #[test]
    fn test_dihedral_angle_trans() {
        let beads = vec![
            make_bead([0.0, 1.0, 0.0]),
            make_bead([0.0, 0.0, 0.0]),
            make_bead([1.0, 0.0, 0.0]),
            make_bead([1.0, -1.0, 0.0]),
        ];
        let dih = CgDihedral::new(0, 1, 2, 3, 0.0, 1.0, 1);
        let phi = dih.dihedral_angle(&beads);
        assert!((phi.abs() - PI).abs() < 1e-10, "Expected pi, got {}", phi);
    }
    #[test]
    fn test_go_model_native_contacts() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let go = GoModel::new(positions, 0.5, 1.0);
        let contacts = go.native_contacts();
        assert_eq!(contacts.len(), 0, "no long-range contacts expected");
    }
    #[test]
    fn test_go_model_native_contact_energy_at_native() {
        let epsilon = 1.0;
        let r0 = 0.4;
        let e = GoModel::native_contact_energy(epsilon, r0, r0);
        assert!(
            (e + epsilon).abs() < 1e-10,
            "energy at r0 = {e}, expected -1.0"
        );
    }
    #[test]
    fn test_go_model_q_value_native() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.38, 0.0, 0.0],
            [0.76, 0.0, 0.0],
            [1.14, 0.0, 0.0],
        ];
        let go = GoModel::new(positions.clone(), 0.9, 1.0);
        let q = go.q_value(&positions);
        assert!((q - 1.0).abs() < 1e-10 || q == 0.0, "Q at native = {q}");
    }
    #[test]
    fn test_go_model_repulsive_energy_decreases_with_distance() {
        let sigma = 0.4;
        let e_close = GoModel::repulsive_energy(sigma, 0.2);
        let e_far = GoModel::repulsive_energy(sigma, 0.8);
        assert!(e_close > e_far, "repulsion should decrease with distance");
    }
    #[test]
    fn test_bead_mapping_centroid() {
        let mapping = BeadMapping::new(vec![vec![0, 1], vec![2, 3]]);
        let aa_positions = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [2.0, 4.0, 0.0],
        ];
        let cg_pos = mapping.map_positions(&aa_positions);
        assert_eq!(cg_pos.len(), 2);
        assert!(
            (cg_pos[0][0] - 1.0).abs() < 1e-12,
            "centroid x = {}",
            cg_pos[0][0]
        );
        assert!(
            (cg_pos[1][1] - 4.0).abs() < 1e-12,
            "centroid y = {}",
            cg_pos[1][1]
        );
    }
    #[test]
    fn test_bead_mapping_masses() {
        let mapping = BeadMapping::new(vec![vec![0, 1], vec![2]]);
        let masses = vec![12.0, 16.0, 1.0];
        let bead_masses = mapping.bead_masses(&masses);
        assert_eq!(bead_masses.len(), 2);
        assert!((bead_masses[0] - 28.0).abs() < 1e-12);
        assert!((bead_masses[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_bead_mapping_velocity_conservation() {
        let mapping = BeadMapping::new(vec![vec![0, 1]]);
        let masses = vec![1.0, 1.0];
        let vels = vec![[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
        let cg_vel = mapping.map_velocities(&vels, &masses);
        assert!(cg_vel[0][0].abs() < 1e-12, "COM vel = {}", cg_vel[0][0]);
    }
    #[test]
    fn test_cg_pair_potential_lj_energy() {
        let sigma = 0.47;
        let epsilon = 5.6;
        let pot =
            CgPairPotential::from_function(0.3, 2.0, 200, |r| MartiniLj::energy(r, sigma, epsilon));
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let e = pot.energy(r_min);
        assert!((e + epsilon).abs() < 0.5, "LJ energy at r_min = {e}");
    }
    #[test]
    fn test_cg_pair_potential_energy_at_endpoints() {
        let pot = CgPairPotential::from_function(0.5, 2.0, 100, |r| r * r);
        let e = pot.energy(10.0);
        assert!(e.is_finite(), "energy at far r = {e}");
    }
    #[test]
    fn test_cg_pair_potential_force_sign() {
        let pot = CgPairPotential::from_function(0.5, 2.0, 200, |r| r * r);
        let f = pot.force(1.0);
        assert!(
            f < 0.0,
            "force should be negative (repulsive) for V=r^2, got {f}"
        );
    }
    #[test]
    fn test_cg_pair_potential_n_points() {
        let pot = CgPairPotential::from_function(0.0, 1.0, 50, |r| r);
        assert_eq!(pot.n_points(), 50);
    }
    #[test]
    fn test_langevin_velocity_scale_less_than_one() {
        let thermo = LangevinThermostat::new(300.0, 1.0);
        let scale = thermo.velocity_scale(0.01);
        assert!(scale < 1.0 && scale > 0.0, "scale = {scale}");
    }
    #[test]
    fn test_langevin_velocity_scale_zero_friction() {
        let thermo = LangevinThermostat::new(300.0, 0.0);
        let scale = thermo.velocity_scale(0.01);
        assert!(
            (scale - 1.0).abs() < 1e-12,
            "zero friction should give scale=1, got {scale}"
        );
    }
    #[test]
    fn test_langevin_random_force_positive_sigma() {
        let thermo = LangevinThermostat::new(300.0, 1.0);
        let sigma = thermo.random_force_sigma(72.0, 0.001);
        assert!(sigma > 0.0, "sigma = {sigma}");
    }
    #[test]
    fn test_cg_cosine_angle_energy_at_equilibrium() {
        let theta0 = PI / 2.0;
        let ka = 50.0;
        let angle = CgCosineAngle::new(0, 1, 2, theta0, ka);
        let beads = vec![
            make_bead([1.0, 0.0, 0.0]),
            make_bead([0.0, 0.0, 0.0]),
            make_bead([0.0, 1.0, 0.0]),
        ];
        let e = angle.energy(&beads);
        assert!(e.abs() < 1e-10, "energy at equilibrium = {e}");
    }
    #[test]
    fn test_cg_cosine_angle_energy_nonzero_displaced() {
        let theta0 = PI / 2.0;
        let ka = 50.0;
        let angle = CgCosineAngle::new(0, 1, 2, theta0, ka);
        let beads = vec![
            make_bead([1.0, 0.0, 0.0]),
            make_bead([0.0, 0.0, 0.0]),
            make_bead([-1.0, 0.0, 0.0]),
        ];
        let e = angle.energy(&beads);
        assert!((e - 25.0).abs() < 1e-8, "energy at 180 deg = {e}");
    }
    #[test]
    fn test_martini_bead_types_count() {
        let types = MartiniBeadType::standard_types();
        assert_eq!(types.len(), 5, "MARTINI standard has 5 bead types");
    }
    #[test]
    fn test_martini_bead_types_mass() {
        for t in MartiniBeadType::standard_types() {
            assert!(
                (t.mass - 72.0).abs() < 1e-10,
                "{} mass = {}",
                t.name,
                t.mass
            );
        }
    }
    #[test]
    fn test_martini_interaction_c_level_rmin() {
        let mi = MartiniInteraction::c_level();
        let r_min = mi.r_min();
        let expected = 2.0_f64.powf(1.0 / 6.0) * 0.47;
        assert!((r_min - expected).abs() < 1e-10, "r_min = {r_min}");
    }
    #[test]
    fn test_martini_interaction_energy_at_rmin() {
        let mi = MartiniInteraction::c_level();
        let e = mi.energy(mi.r_min());
        assert!(
            (e + mi.epsilon).abs() < 0.01,
            "energy at r_min should be -epsilon: {e}"
        );
    }
    #[test]
    fn test_martini_interaction_force_repulsive_close() {
        let mi = MartiniInteraction::c_level();
        let f = mi.force_magnitude(0.1);
        assert!(f > 0.0, "repulsive force at close range: {f}");
    }
    #[test]
    fn test_martini_interaction_d_level_weaker_than_c() {
        let c = MartiniInteraction::c_level();
        let d = MartiniInteraction::d_level();
        let r = c.r_min();
        assert!(
            d.energy(r) > c.energy(r),
            "D level weaker (higher energy at r_min)"
        );
    }
    #[test]
    fn test_shifted_martini_lj_zero_at_cutoff() {
        let mi = MartiniInteraction::c_level();
        let r_cut = 1.2;
        let slj = ShiftedMartiniLj::new(mi, r_cut);
        let e = slj.energy(r_cut);
        assert!(e.abs() < 1e-10, "shifted LJ at r_cut should be 0, got {e}");
    }
    #[test]
    fn test_shifted_martini_lj_zero_beyond_cutoff() {
        let mi = MartiniInteraction::c_level();
        let slj = ShiftedMartiniLj::new(mi, 1.2);
        assert_eq!(slj.energy(2.0), 0.0);
    }
    #[test]
    fn test_enm_n_contacts_two_close() {
        let nodes = vec![
            EnmNode::new([0.0, 0.0, 0.0], "A"),
            EnmNode::new([0.5, 0.0, 0.0], "B"),
        ];
        let enm = EnmNodeModel::new(nodes, 100.0, 1.0);
        assert_eq!(enm.n_contacts(), 1);
    }
    #[test]
    fn test_enm_n_contacts_far() {
        let nodes = vec![
            EnmNode::new([0.0, 0.0, 0.0], "A"),
            EnmNode::new([5.0, 0.0, 0.0], "B"),
        ];
        let enm = EnmNodeModel::new(nodes, 100.0, 1.0);
        assert_eq!(enm.n_contacts(), 0);
    }
    #[test]
    fn test_enm_displaced_energy_at_native() {
        let ref_nodes = vec![
            EnmNode::new([0.0, 0.0, 0.0], "A"),
            EnmNode::new([0.5, 0.0, 0.0], "B"),
        ];
        let enm = EnmNodeModel::new(ref_nodes.clone(), 100.0, 1.0);
        let e = enm.displaced_energy(&ref_nodes);
        assert!(e.abs() < 1e-12, "displaced energy at native = {e}");
    }
    #[test]
    fn test_enm_displaced_energy_increases() {
        let ref_nodes = vec![
            EnmNode::new([0.0, 0.0, 0.0], "A"),
            EnmNode::new([0.5, 0.0, 0.0], "B"),
        ];
        let displaced_nodes = vec![
            EnmNode::new([0.0, 0.0, 0.0], "A"),
            EnmNode::new([0.8, 0.0, 0.0], "B"),
        ];
        let enm = EnmNodeModel::new(displaced_nodes, 100.0, 1.0);
        let e = enm.displaced_energy(&ref_nodes);
        assert!(e > 0.0, "displaced energy should be positive: {e}");
    }
    #[test]
    fn test_enm_connectivity() {
        let nodes = vec![
            EnmNode::new([0.0, 0.0, 0.0], "A"),
            EnmNode::new([0.5, 0.0, 0.0], "B"),
            EnmNode::new([1.0, 0.0, 0.0], "C"),
        ];
        let enm = EnmNodeModel::new(nodes, 100.0, 0.7);
        let conn = enm.connectivity();
        assert_eq!(conn[0], 1);
        assert_eq!(conn[1], 2);
        assert_eq!(conn[2], 1);
    }
    #[test]
    fn test_cg_water_bead_mass() {
        assert!((CgWaterBead::mass() - 72.06).abs() < 1e-6);
    }
    #[test]
    fn test_cg_water_energy_at_rmin() {
        let r_min = 2.0_f64.powf(1.0 / 6.0) * 0.47;
        let e = cg_water_energy(r_min);
        assert!((e + 5.0).abs() < 0.01, "energy at r_min = {e}");
    }
    #[test]
    fn test_build_cg_water_box_count() {
        let beads = build_cg_water_box(3, 0.5);
        assert_eq!(beads.len(), 27, "3^3 = 27 beads");
    }
    #[test]
    fn test_build_cg_water_box_positions() {
        let beads = build_cg_water_box(2, 0.47);
        assert!((beads[0].pos[0]).abs() < 1e-12);
        assert!((beads.last().unwrap().pos[0] - 0.47).abs() < 1e-10);
    }
    #[test]
    fn test_backmap_bead_single_atom() {
        let pos = [1.0, 2.0, 3.0];
        let atoms = backmap_bead(pos, 1, 0.1);
        assert_eq!(atoms.len(), 1);
        assert!((atoms[0][0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_backmap_bead_four_atoms() {
        let pos = [0.0, 0.0, 0.0];
        let atoms = backmap_bead(pos, 4, 0.1);
        assert_eq!(atoms.len(), 4);
        for a in &atoms {
            let d = length(sub(*a, pos));
            assert!(d <= 0.11, "atom at distance {d} from bead");
        }
    }
    #[test]
    fn test_backmap_molecule_count() {
        let beads = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let atoms = backmap_molecule(&beads, 4, 0.1);
        assert_eq!(atoms.len(), 8, "2 beads × 4 atoms = 8");
    }
    #[test]
    fn test_backmap_bead_empty() {
        let pos = [0.0, 0.0, 0.0];
        let atoms = backmap_bead(pos, 0, 0.1);
        assert!(atoms.is_empty());
    }
    #[test]
    fn test_cg_fourier_dihedral_angle_zero() {
        let beads = vec![
            make_bead([0.0, 1.0, 0.0]),
            make_bead([0.0, 0.0, 0.0]),
            make_bead([1.0, 0.0, 0.0]),
            make_bead([1.0, 1.0, 0.0]),
        ];
        let phi = CgFourierDihedral::dihedral_angle(&beads, 0, 1, 2, 3);
        assert!(phi.abs() < 1e-10, "coplanar dihedral = {phi}");
    }
    #[test]
    fn test_cg_fourier_dihedral_energy_at_zero() {
        let beads = vec![
            make_bead([0.0, 1.0, 0.0]),
            make_bead([0.0, 0.0, 0.0]),
            make_bead([1.0, 0.0, 0.0]),
            make_bead([1.0, 1.0, 0.0]),
        ];
        let dih = CgFourierDihedral::new(0, 1, 2, 3, vec![1.0], vec![0.0]);
        let e = dih.energy(&beads);
        assert!((e - 2.0).abs() < 1e-8, "dihedral energy = {e}");
    }
    #[test]
    fn test_cg_rg_single_bead() {
        let beads = vec![make_bead([1.0, 2.0, 3.0])];
        let rg = cg_radius_of_gyration(&beads);
        assert!(rg.abs() < 1e-12, "Rg of single bead = {rg}");
    }
    #[test]
    fn test_cg_rg_two_equal_beads() {
        let beads = vec![make_bead([-1.0, 0.0, 0.0]), make_bead([1.0, 0.0, 0.0])];
        let rg = cg_radius_of_gyration(&beads);
        assert!((rg - 1.0).abs() < 1e-10, "Rg = {rg}, expected 1.0");
    }
    #[test]
    fn test_cg_rg_empty() {
        let beads: Vec<CgBead> = vec![];
        assert_eq!(cg_radius_of_gyration(&beads), 0.0);
    }
    #[test]
    fn test_end_to_end_distance() {
        let beads = vec![make_bead([0.0, 0.0, 0.0]), make_bead([3.0, 4.0, 0.0])];
        let d = end_to_end_distance(&beads);
        assert!((d - 5.0).abs() < 1e-10, "end-to-end = {d}");
    }
    #[test]
    fn test_end_to_end_single() {
        let beads = vec![make_bead([0.0, 0.0, 0.0])];
        let d = end_to_end_distance(&beads);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_cg_kinetic_energy() {
        let mut b = make_bead([0.0, 0.0, 0.0]);
        b.vel = [1.0, 0.0, 0.0];
        let ke = cg_kinetic_energy(&[b]);
        assert!((ke - 0.5 * 72.0).abs() < 1e-10, "KE = {ke}");
    }
    #[test]
    fn test_cg_kinetic_energy_zero() {
        let b = make_bead([0.0, 0.0, 0.0]);
        let ke = cg_kinetic_energy(&[b]);
        assert_eq!(ke, 0.0);
    }
}
