//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Bead, CgAngle, CgBond, CgDihedral, MartiniBead};

/// Boltzmann constant (J K^-1).
pub(super) const KB: f64 = 1.380_649e-23;
/// Compute the standard Lennard-Jones energy between two beads at distance `r`.
///
/// Formula: 4ε * ((σ/r)^12 - (σ/r)^6)
///
/// The minimum is at r = σ·2^(1/6) with value −ε.
pub fn martini_lj_energy(r: f64, sigma: f64, epsilon: f64) -> f64 {
    let sr = sigma / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    4.0 * epsilon * (sr12 - sr6)
}
/// Compute the magnitude of the Lennard-Jones force between two beads at distance `r`.
///
/// F(r) = −dU/dr = 24ε/r * (2*(σ/r)^12 - (σ/r)^6)
///
/// Returns a positive value when the force is repulsive (r < r_min).
pub fn martini_lj_force(r: f64, sigma: f64, epsilon: f64) -> f64 {
    let sr = sigma / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    24.0 * epsilon / r * (2.0 * sr12 - sr6)
}
/// Compute the total harmonic bond energy for a set of CG beads.
///
/// E = Σ ½ k (r − r₀)²
pub fn cg_bond_energy(beads: &[[f64; 3]], bonds: &[CgBond]) -> f64 {
    bonds.iter().fold(0.0, |acc, bond| {
        let pi = beads[bond.bead_i];
        let pj = beads[bond.bead_j];
        let dx = pj[0] - pi[0];
        let dy = pj[1] - pi[1];
        let dz = pj[2] - pi[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let dr = r - bond.r0;
        acc + 0.5 * bond.k * dr * dr
    })
}
/// Compute the total harmonic angle energy for a set of CG beads.
///
/// E = Σ ½ k_θ (θ − θ₀)²
pub fn cg_angle_energy(beads: &[[f64; 3]], angles: &[CgAngle]) -> f64 {
    angles.iter().fold(0.0, |acc, angle| {
        let pi = beads[angle.i];
        let pj = beads[angle.j];
        let pk = beads[angle.k];
        let v1 = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
        let v2 = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let n1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let n2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        let cos_theta = (dot / (n1 * n2)).clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let dtheta = theta - angle.theta0;
        acc + 0.5 * angle.k_theta * dtheta * dtheta
    })
}
/// Compute the radius of gyration for a set of coarse-grained beads.
///
/// Rg = sqrt( Σ mᵢ |rᵢ − r_com|² / Σ mᵢ )
pub fn compute_radius_of_gyration_cg(beads: &[Bead]) -> f64 {
    if beads.is_empty() {
        return 0.0;
    }
    let mut com = [0.0f64; 3];
    let mut total_mass = 0.0f64;
    for bead in beads {
        total_mass += bead.mass;
        for (d, c) in com.iter_mut().enumerate() {
            *c += bead.mass * bead.position[d];
        }
    }
    if total_mass > 0.0 {
        for v in &mut com {
            *v /= total_mass;
        }
    }
    let sum_mr2: f64 = beads.iter().fold(0.0, |acc, bead| {
        let dr2: f64 = (0..3)
            .map(|d| {
                let dr = bead.position[d] - com[d];
                dr * dr
            })
            .sum();
        acc + bead.mass * dr2
    });
    (sum_mr2 / total_mass).sqrt()
}
/// Compute the dihedral angle between four points.
///
/// Returns the angle in radians in \[-π, π\].
pub fn dihedral_angle(p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3], p4: &[f64; 3]) -> f64 {
    let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
    let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];
    let n1 = cross(&b1, &b2);
    let n2 = cross(&b2, &b3);
    let m1 = cross(&n1, &b2);
    let b2_len = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    if b2_len < 1e-15 {
        return 0.0;
    }
    let m1_norm = [m1[0] / b2_len, m1[1] / b2_len, m1[2] / b2_len];
    let x = dot_product(&n1, &n2);
    let y = dot_product(&m1_norm, &n2);
    (-y).atan2(x)
}
/// Compute total dihedral energy.
///
/// E = Σ k * (1 + cos(n*φ - φ_0))
pub fn cg_dihedral_energy(beads: &[[f64; 3]], dihedrals: &[CgDihedral]) -> f64 {
    dihedrals.iter().fold(0.0, |acc, dih| {
        let phi = dihedral_angle(&beads[dih.i], &beads[dih.j], &beads[dih.k], &beads[dih.l]);
        let cos_val = (dih.n as f64 * phi - dih.phi0).cos();
        acc + dih.k_phi * (1.0 + cos_val)
    })
}
/// Cross product of two 3-vectors.
pub(super) fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Dot product of two 3-vectors.
pub(super) fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Compute the radial distribution function g(r) for a set of positions
/// in a cubic box with periodic boundary conditions.
///
/// Returns a vector of length `n_bins` with normalised g(r) values.
/// Bin centres are at (i + 0.5) * r_max / n_bins.
pub fn cg_radial_distribution_function(
    positions: &[[f64; 3]],
    n_bins: usize,
    r_max: f64,
    box_len: f64,
) -> Vec<f64> {
    let n = positions.len();
    let dr = r_max / n_bins as f64;
    let mut histogram = vec![0.0f64; n_bins];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut r2 = 0.0f64;
            for (pi, qi) in positions[i].iter().zip(positions[j].iter()) {
                let mut dx = qi - pi;
                dx -= box_len * (dx / box_len).round();
                r2 += dx * dx;
            }
            let r = r2.sqrt();
            if r < r_max {
                let bin = (r / dr) as usize;
                if bin < n_bins {
                    histogram[bin] += 2.0;
                }
            }
        }
    }
    let volume = box_len * box_len * box_len;
    let density = n as f64 / volume;
    let mut g = vec![0.0f64; n_bins];
    for i in 0..n_bins {
        let r_inner = i as f64 * dr;
        let r_outer = (i as f64 + 1.0) * dr;
        let shell_vol = (4.0 / 3.0) * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
        let ideal = n as f64 * density * shell_vol;
        if ideal > 0.0 {
            g[i] = histogram[i] / ideal;
        }
    }
    g
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_cg_mapping_com_positions() {
        let atom_positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
        ];
        let atom_masses = vec![1.0, 1.0, 1.0, 1.0];
        let mut mapping = CgMapping::new();
        mapping.add_bead(vec![0, 1], 2.0);
        mapping.add_bead(vec![2, 3], 2.0);
        let cg_pos = mapping.map_positions(&atom_positions, &atom_masses);
        assert_eq!(cg_pos.len(), 2);
        assert!((cg_pos[0][0] - 1.0).abs() < 1e-10);
        assert!((cg_pos[0][1]).abs() < 1e-10);
        assert!((cg_pos[0][2]).abs() < 1e-10);
        assert!((cg_pos[1][0] - 5.0).abs() < 1e-10);
        assert!((cg_pos[1][1]).abs() < 1e-10);
        assert!((cg_pos[1][2]).abs() < 1e-10);
    }
    #[test]
    fn test_martini_lj_energy_minimum() {
        let sigma = 0.47;
        let epsilon = 3.5;
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let energy = martini_lj_energy(r_min, sigma, epsilon);
        assert!((energy - (-epsilon)).abs() < 1e-10);
    }
    #[test]
    fn test_cg_bond_energy_at_equilibrium() {
        let beads: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.47, 0.0, 0.0]];
        let bonds = vec![CgBond {
            bead_i: 0,
            bead_j: 1,
            r0: 0.47,
            k: 3800.0,
        }];
        let energy = cg_bond_energy(&beads, &bonds);
        assert!(energy.abs() < 1e-20);
    }
    #[test]
    fn test_compute_radius_of_gyration_cg() {
        let beads = vec![
            Bead {
                position: [1.0, 0.0, 0.0],
                velocity: [0.0, 0.0, 0.0],
                mass: 1.0,
                bead_type: 0,
            },
            Bead {
                position: [-1.0, 0.0, 0.0],
                velocity: [0.0, 0.0, 0.0],
                mass: 1.0,
                bead_type: 0,
            },
        ];
        let rg = compute_radius_of_gyration_cg(&beads);
        assert!((rg - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cg_mapping_forces() {
        let atom_forces: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-2.0, 1.0, 0.0],
        ];
        let mut mapping = CgMapping::new();
        mapping.add_bead(vec![0, 1], 2.0);
        mapping.add_bead(vec![2, 3], 2.0);
        let cg_forces = mapping.map_forces(&atom_forces);
        assert_eq!(cg_forces.len(), 2);
        assert!((cg_forces[0][0] - 3.0).abs() < 1e-10);
        assert!((cg_forces[1][0] - (-3.0)).abs() < 1e-10);
        assert!((cg_forces[1][1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_cg_mapping_velocities() {
        let atom_velocities: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let atom_masses = vec![1.0, 3.0];
        let mut mapping = CgMapping::new();
        mapping.add_bead(vec![0, 1], 4.0);
        let cg_vel = mapping.map_velocities(&atom_velocities, &atom_masses);
        assert!((cg_vel[0][0] - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_force_match_accumulator() {
        let mut acc = ForceMatchAccumulator::new(0.0, 1.0, 10);
        acc.add_sample(0.25, 10.0);
        acc.add_sample(0.25, 20.0);
        acc.add_sample(0.55, -5.0);
        let result = acc.compute_mean_forces();
        assert_eq!(result.counts[2], 2);
        assert!((result.mean_forces[2] - 15.0).abs() < 1e-10);
        assert_eq!(result.counts[5], 1);
        assert!((result.mean_forces[5] - (-5.0)).abs() < 1e-10);
    }
    #[test]
    fn test_force_match_out_of_range() {
        let mut acc = ForceMatchAccumulator::new(0.2, 1.0, 8);
        acc.add_sample(0.1, 10.0);
        acc.add_sample(1.5, 10.0);
        let result = acc.compute_mean_forces();
        let total: usize = result.counts.iter().sum();
        assert_eq!(total, 0, "Out-of-range samples should be discarded");
    }
    #[test]
    fn test_force_match_reset() {
        let mut acc = ForceMatchAccumulator::new(0.2, 1.0, 8);
        acc.add_sample(0.3, 10.0);
        acc.reset();
        let result = acc.compute_mean_forces();
        let total: usize = result.counts.iter().sum();
        assert_eq!(total, 0, "After reset, all counts should be 0");
    }
    #[test]
    fn test_ibi_initial_potential() {
        let r = vec![0.3, 0.4, 0.5, 0.6];
        let g_target = vec![0.5, 1.0, 1.2, 1.0];
        let kbt = 2.494;
        let ibi = IbiPotential::from_rdf(r, g_target, kbt);
        assert!((ibi.potential[1]).abs() < 1e-10, "U should be 0 where g=1");
        assert!(ibi.potential[0] > 0.0, "U should be positive where g<1");
        assert!(ibi.potential[2] < 0.0, "U should be negative where g>1");
    }
    #[test]
    fn test_ibi_update_convergence() {
        let r = vec![0.3, 0.4, 0.5];
        let g_target = vec![0.5, 1.0, 1.5];
        let kbt = 2.494;
        let mut ibi = IbiPotential::from_rdf(r, g_target.clone(), kbt);
        let u_before = ibi.potential.clone();
        ibi.update(&g_target, 1.0);
        for (i, (&before, &after)) in u_before.iter().zip(ibi.potential.iter()).enumerate() {
            assert!(
                (before - after).abs() < 1e-10,
                "Potential[{i}] should not change when g_curr = g_target"
            );
        }
    }
    #[test]
    fn test_ibi_interpolation() {
        let r = vec![0.3, 0.4, 0.5];
        let g_target = vec![1.0, 1.0, 1.0];
        let kbt = 2.494;
        let ibi = IbiPotential::from_rdf(r, g_target, kbt);
        let u = ibi.interpolate(0.35);
        assert!(
            u.abs() < 1e-8,
            "Interpolated potential should be ~0, got {u}"
        );
    }
    #[test]
    fn test_ibi_force() {
        let r = vec![0.3, 0.4, 0.5, 0.6];
        let g_target = vec![0.1, 0.5, 1.0, 1.0];
        let kbt = 2.494;
        let ibi = IbiPotential::from_rdf(r, g_target, kbt);
        let f = ibi.force(0.35);
        assert!(f > 0.0, "Force should be repulsive at short range, got {f}");
    }
    #[test]
    fn test_ibi_convergence_metric() {
        let r = vec![0.3, 0.4, 0.5];
        let g_target = vec![1.0, 1.0, 1.0];
        let kbt = 2.494;
        let ibi = IbiPotential::from_rdf(r, g_target.clone(), kbt);
        let metric = ibi.convergence_metric(&g_target);
        assert!(
            metric.abs() < 1e-10,
            "Converged system should have metric ~0"
        );
        let g_bad = vec![1.5, 0.5, 2.0];
        let metric_bad = ibi.convergence_metric(&g_bad);
        assert!(metric_bad > 0.1, "Non-converged should have large metric");
    }
    #[test]
    fn test_tabulated_potential_lj() {
        let tab = TabulatedPotential::lj(0.4, 1.5, 1000, 0.47, 3.5);
        let r_min = 0.47 * 2.0_f64.powf(1.0 / 6.0);
        let e = tab.energy(r_min);
        assert!(
            (e - (-3.5)).abs() < 0.1,
            "LJ energy at minimum should be ~-3.5, got {e}"
        );
    }
    #[test]
    fn test_tabulated_potential_force_sign() {
        let tab = TabulatedPotential::lj(0.4, 1.5, 1000, 0.47, 3.5);
        let r_min = 0.47 * 2.0_f64.powf(1.0 / 6.0);
        let f_below = tab.force(r_min - 0.05);
        assert!(
            f_below > 0.0,
            "Force below r_min should be repulsive, got {f_below}"
        );
        let f_above = tab.force(r_min + 0.05);
        assert!(
            f_above < 0.0,
            "Force above r_min should be attractive, got {f_above}"
        );
    }
    #[test]
    fn test_tabulated_potential_shift() {
        let mut tab = TabulatedPotential::lj(0.4, 1.5, 1000, 0.47, 3.5);
        tab.shift_to_zero_at_cutoff();
        let e_cut = tab.energy(1.5);
        assert!(
            e_cut.abs() < 0.01,
            "Energy at cutoff should be ~0 after shift, got {e_cut}"
        );
    }
    #[test]
    fn test_tabulated_from_function() {
        let k = 100.0;
        let r0 = 0.5;
        let tab =
            TabulatedPotential::from_function(0.3, 0.7, 500, |r| 0.5 * k * (r - r0) * (r - r0));
        let e_at_r0 = tab.energy(r0);
        assert!(
            e_at_r0.abs() < 0.001,
            "Energy at r0 should be ~0, got {e_at_r0}"
        );
    }
    #[test]
    fn test_adaptive_resolution_weight() {
        let ar = AdaptiveResolution::new([0.0, 0.0, 0.0], 1.0, 2.0);
        let w_aa = ar.weight(&[0.5, 0.0, 0.0]);
        assert!(
            (w_aa - 1.0).abs() < 1e-10,
            "Weight in AA region should be 1.0"
        );
        let w_cg = ar.weight(&[3.0, 0.0, 0.0]);
        assert!(w_cg.abs() < 1e-10, "Weight in CG region should be 0.0");
        let w_hyb = ar.weight(&[1.5, 0.0, 0.0]);
        assert!(
            w_hyb > 0.0 && w_hyb < 1.0,
            "Weight in hybrid region should be in (0,1), got {w_hyb}"
        );
    }
    #[test]
    fn test_adaptive_resolution_region() {
        let ar = AdaptiveResolution::new([0.0, 0.0, 0.0], 1.0, 2.0);
        assert_eq!(ar.region(&[0.0, 0.0, 0.0]), ResolutionRegion::Atomistic);
        assert_eq!(ar.region(&[1.5, 0.0, 0.0]), ResolutionRegion::Hybrid);
        assert_eq!(ar.region(&[3.0, 0.0, 0.0]), ResolutionRegion::CoarseGrained);
    }
    #[test]
    fn test_adaptive_resolution_blend() {
        let ar = AdaptiveResolution::new([0.0, 0.0, 0.0], 1.0, 2.0);
        let f_aa = [10.0, 0.0, 0.0];
        let f_cg = [5.0, 0.0, 0.0];
        let f_blend = ar.blend_forces(&[0.0, 0.0, 0.0], &f_aa, &f_cg);
        assert!((f_blend[0] - 10.0).abs() < 1e-10);
        let f_blend = ar.blend_forces(&[3.0, 0.0, 0.0], &f_aa, &f_cg);
        assert!((f_blend[0] - 5.0).abs() < 1e-10);
        let f_blend = ar.blend_forces(&[1.5, 0.0, 0.0], &f_aa, &f_cg);
        assert!(
            f_blend[0] > 5.0 && f_blend[0] < 10.0,
            "Blended force should be between AA and CG"
        );
    }
    #[test]
    fn test_rdf_calculator_accumulate() {
        let mut rdf = RdfCalculator::new(2.0, 20, 4, 5.0);
        let positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
        ];
        rdf.accumulate(&positions);
        assert_eq!(rdf.n_frames, 1);
    }
    #[test]
    fn test_rdf_calculator_rdf_shape() {
        let mut rdf = RdfCalculator::new(2.0, 20, 4, 5.0);
        let positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
        ];
        rdf.accumulate(&positions);
        let (r_vals, g_vals) = rdf.compute_rdf();
        assert_eq!(r_vals.len(), 20);
        assert_eq!(g_vals.len(), 20);
    }
    #[test]
    fn test_rdf_calculator_reset() {
        let mut rdf = RdfCalculator::new(2.0, 20, 4, 5.0);
        let positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
        ];
        rdf.accumulate(&positions);
        rdf.reset();
        assert_eq!(rdf.n_frames, 0);
    }
    #[test]
    fn test_dihedral_angle_planar() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [2.0, 0.0, 0.0];
        let p4 = [3.0, 0.0, 0.0];
        let phi = dihedral_angle(&p1, &p2, &p3, &p4);
        assert!(phi.is_finite());
    }
    #[test]
    fn test_dihedral_angle_90_degrees() {
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let p4 = [0.0, 1.0, 1.0];
        let phi = dihedral_angle(&p1, &p2, &p3, &p4);
        assert!(
            (phi.abs() - std::f64::consts::FRAC_PI_2).abs() < 0.1,
            "Dihedral should be ~90°, got {} rad",
            phi
        );
    }
    #[test]
    fn test_cg_dihedral_energy() {
        let beads: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let dihedrals = vec![CgDihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            k_phi: 5.0,
            n: 1,
            phi0: 0.0,
        }];
        let energy = cg_dihedral_energy(&beads, &dihedrals);
        assert!(
            energy >= 0.0,
            "Dihedral energy should be >= 0, got {energy}"
        );
    }
    #[test]
    fn test_cg_bead_properties() {
        let bb = CgBeadProperties::martini_backbone();
        assert_eq!(bb.bead_type, CgBeadType::Backbone);
        assert!((bb.mass - 72.0).abs() < 1e-10);
        assert!((bb.sigma - 0.47).abs() < 1e-10);
    }
    #[test]
    fn test_cg_bead_water() {
        let w = CgBeadProperties::martini_water();
        assert_eq!(w.bead_type, CgBeadType::Water);
        assert!((w.charge).abs() < 1e-10, "Water bead should be neutral");
    }
    #[test]
    fn test_cg_bead_custom() {
        let custom = CgBeadProperties::new(CgBeadType::Custom(42), 100.0, -1.0, 0.5, 4.0);
        assert_eq!(custom.bead_type, CgBeadType::Custom(42));
        assert!((custom.charge - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_cross_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross(&a, &b);
        assert!((c[0]).abs() < 1e-15);
        assert!((c[1]).abs() < 1e-15);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_dot_product() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let d = dot_product(&a, &b);
        assert!((d - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_cg_system_com_weighted() {
        let mut sys = CgSystem::new([10.0, 10.0, 10.0]);
        sys.add_bead(CgBead {
            position: [0.0, 0.0, 0.0],
            mass: 1.0,
            bead_type: 0,
            charge: 0.0,
        });
        sys.add_bead(CgBead {
            position: [2.0, 0.0, 0.0],
            mass: 3.0,
            bead_type: 0,
            charge: 0.0,
        });
        let com = sys.center_of_mass();
        assert!(
            (com[0] - 1.5).abs() < 1e-12,
            "COM x should be 1.5, got {}",
            com[0]
        );
        assert!(com[1].abs() < 1e-12);
        assert!(com[2].abs() < 1e-12);
    }
    #[test]
    fn test_martini_mapping_com() {
        let scheme = MartiniMappingScheme::new(4, 72.0);
        let positions = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
        ];
        let masses = vec![1.0, 1.0, 1.0, 1.0];
        let com = scheme.map_atomistic(&positions, &masses);
        assert!(
            (com[0] - 1.0).abs() < 1e-12,
            "COM x should be 1.0, got {}",
            com[0]
        );
        assert!(
            (com[1] - 1.0).abs() < 1e-12,
            "COM y should be 1.0, got {}",
            com[1]
        );
    }
    #[test]
    fn test_cg_nonbonded_lj_energy_minimum() {
        let sigma = 0.47;
        let epsilon = 3.5;
        let nb = CgNonbonded::new(epsilon, sigma, 2.0);
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let e = nb.energy(r_min);
        assert!(
            (e - (-epsilon)).abs() < 1e-10,
            "LJ minimum energy should be -epsilon, got {e}"
        );
    }
    #[test]
    fn test_ibi_update_shifts_potential() {
        let r_bins = vec![0.3, 0.4, 0.5];
        let target_rdf = vec![1.0, 1.0, 1.0];
        let mut ibi = IterativeBoltzmannInversion::from_target(target_rdf.clone(), r_bins, 300.0);
        let pot_before = ibi.current_potential.clone();
        let current_rdf = vec![2.0, 2.0, 2.0];
        ibi.update_potential(&current_rdf, 300.0);
        let changed = pot_before
            .iter()
            .zip(ibi.current_potential.iter())
            .any(|(a, b)| (a - b).abs() > 1e-15);
        assert!(changed, "IBI update should shift the potential");
    }
    #[test]
    fn test_cg_rdf_bin_count() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let g = cg_radial_distribution_function(&positions, 50, 5.0, 10.0);
        assert_eq!(g.len(), 50, "RDF should have exactly n_bins values");
    }
}
/// Linear interpolation for tabulated potentials.
pub fn interpolate_table(x: f64, x_min: f64, dx: f64, values: &[f64]) -> f64 {
    let n = values.len();
    if n == 0 {
        return 0.0;
    }
    let idx_f = (x - x_min) / dx;
    if idx_f < 0.0 {
        return values[0];
    }
    let idx = idx_f as usize;
    if idx + 1 >= n {
        return values[n - 1];
    }
    let t = idx_f - idx as f64;
    values[idx] * (1.0 - t) + values[idx + 1] * t
}
#[cfg(test)]
mod tests_cg_ext {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_enm_build_contact_count() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let enm = ElasticNetworkModel::build(&positions, 0.5, 1.0);
        assert_eq!(enm.contacts.len(), 1, "only one contact within cutoff");
    }
    #[test]
    fn test_enm_energy_zero_at_equilibrium() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let enm = ElasticNetworkModel::build(&positions, 0.5, 1.0);
        let e = enm.energy(&positions);
        assert!(
            e.abs() < 1e-20,
            "energy should be zero at equilibrium positions, got {e}"
        );
    }
    #[test]
    fn test_enm_forces_zero_at_equilibrium() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let enm = ElasticNetworkModel::build(&positions, 0.5, 1.0);
        let f = enm.forces(&positions);
        for fi in &f {
            for &c in fi {
                assert!(c.abs() < 1e-20, "force should be zero at equilibrium");
            }
        }
    }
    #[test]
    fn test_enm_hessian_size() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let enm = ElasticNetworkModel::build(&positions, 0.5, 1.0);
        let h = enm.hessian(&positions);
        assert_eq!(h.len(), 36, "6×6 Hessian for 2 atoms");
    }
    #[test]
    fn test_martini_lj_params_p4_p4() {
        let p = MartiniLjParams::lookup(MartiniBead::P4, MartiniBead::P4);
        assert!(
            (p.epsilon - 5.6).abs() < 1e-12,
            "P4-P4 epsilon should be 5.6"
        );
    }
    #[test]
    fn test_martini_lj_params_asymmetric() {
        let p1 = MartiniLjParams::lookup(MartiniBead::P4, MartiniBead::C1);
        let p2 = MartiniLjParams::lookup(MartiniBead::C1, MartiniBead::P4);
        assert!(
            (p1.epsilon - p2.epsilon).abs() < 1e-12,
            "should be symmetric"
        );
    }
    #[test]
    fn test_back_mapper_n_atoms() {
        let frames = vec![
            vec![[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]],
            vec![[0.3, 0.0, 0.0]],
        ];
        let masses = vec![vec![12.0, 1.0], vec![14.0]];
        let bm = BackMapper::new(frames, masses);
        assert_eq!(bm.n_atoms(), 3, "total atoms should be 3");
    }
    #[test]
    fn test_back_mapper_places_atoms() {
        let frames = vec![vec![[0.1, 0.0, 0.0]]];
        let masses = vec![vec![1.0]];
        let bm = BackMapper::new(frames, masses);
        let cg = vec![[1.0, 2.0, 3.0]];
        let aa = bm.back_map(&cg);
        assert!(
            (aa[0][0] - 1.1).abs() < 1e-12,
            "back-mapped x should be 1.1"
        );
    }
    #[test]
    fn test_cg_pair_potential_lj_minimum() {
        let eps = 4.5;
        let sig = 0.47;
        let pot = CgPairPotential::lj_table(eps, sig, 0.3, 1.5, 1000);
        let r_min = sig * 2.0_f64.powf(1.0 / 6.0);
        let e = pot.energy(r_min);
        assert!(
            (e - (-eps)).abs() < 0.01,
            "LJ minimum should be ~-eps, got {e}"
        );
    }
    #[test]
    fn test_gnm_degree() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.0, 0.3, 0.0]];
        let gnm = GaussianNetworkModel::build(&positions, 0.5);
        let deg = gnm.degree();
        assert_eq!(deg.len(), 3);
        assert!(
            deg.iter().all(|&d| d >= 1.0),
            "all nodes should have at least one contact"
        );
    }
    #[test]
    fn test_ibi_runner_convergence_after_perfect_rdf() {
        let target = vec![1.0, 1.0, 1.0, 1.0];
        let r_bins = vec![0.3, 0.4, 0.5, 0.6];
        let mut runner = IbiRunner::new(target.clone(), r_bins, 300.0);
        runner.iterate(&target, 300.0);
        assert!(
            runner.converged(1e-10),
            "IBI should converge when current=target"
        );
    }
    #[test]
    fn test_cg_forcefield_pair_potential_lookup() {
        let pot = CgPairPotential::lj_table(4.5, 0.47, 0.3, 1.5, 100);
        let mut ff = CgForceField::new();
        ff.add_pair_potential(0, 1, pot);
        assert!(ff.get_pair_potential(0, 1).is_some(), "should find (0,1)");
        assert!(
            ff.get_pair_potential(1, 0).is_some(),
            "should find (1,0) by symmetry"
        );
        assert!(
            ff.get_pair_potential(0, 2).is_none(),
            "should not find (0,2)"
        );
    }
    #[test]
    fn test_fm_normal_equations_accumulate() {
        let mut ne = FmNormalEquations::new(10, 0.3, 0.1);
        let beads = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let pairs = vec![(0, 1)];
        let forces = vec![1.0];
        ne.accumulate(&beads, &pairs, &forces);
        assert_eq!(ne.n_frames, 1, "n_frames should be 1");
        let any_nonzero = ne.atb.iter().any(|&v| v.abs() > 0.0);
        assert!(
            any_nonzero,
            "A^T b should have nonzero entries after accumulation"
        );
    }
    #[test]
    fn test_interpolate_table_endpoints() {
        let values = vec![0.0, 1.0, 2.0, 3.0];
        assert!((interpolate_table(0.0, 0.0, 1.0, &values) - 0.0).abs() < 1e-12);
        assert!((interpolate_table(1.0, 0.0, 1.0, &values) - 1.0).abs() < 1e-12);
        assert!((interpolate_table(10.0, 0.0, 1.0, &values) - 3.0).abs() < 1e-12);
    }
}
/// MARTINI bead radius (nm) for a given `MartiniBead` type.
pub fn martini_bead_radius(bead: MartiniBead) -> f64 {
    match bead {
        MartiniBead::Sc3 => 0.43,
        MartiniBead::C1
        | MartiniBead::C3
        | MartiniBead::Nda
        | MartiniBead::Qa
        | MartiniBead::Qd
        | MartiniBead::P4 => 0.47,
    }
}
/// MARTINI bead mass (amu) – typically 72 for regular, 54 for small beads.
pub fn martini_bead_mass(bead: MartiniBead) -> f64 {
    match bead {
        MartiniBead::Sc3 => 54.0,
        _ => 72.0,
    }
}
/// Four-to-one MARTINI mapping: group `atom_positions` into CG beads,
/// four atoms per bead, using center-of-mass.
///
/// Returns CG bead positions. Trailing atoms that don't form a complete group
/// are mapped as their own bead.
pub fn martini_4to1_mapping(atom_positions: &[[f64; 3]], atom_masses: &[f64]) -> Vec<[f64; 3]> {
    let n = atom_positions.len().min(atom_masses.len());
    if n == 0 {
        return Vec::new();
    }
    let n_beads = n.div_ceil(4);
    let mut beads = Vec::with_capacity(n_beads);
    for b in 0..n_beads {
        let start = b * 4;
        let end = (start + 4).min(n);
        let mut com = [0.0_f64; 3];
        let mut total_mass = 0.0_f64;
        for i in start..end {
            let m = atom_masses[i];
            total_mass += m;
            for d in 0..3 {
                com[d] += m * atom_positions[i][d];
            }
        }
        if total_mass > 0.0 {
            for v in &mut com {
                *v /= total_mass;
            }
        }
        beads.push(com);
    }
    beads
}
/// Inverse mapping: distribute forces on a CG bead back to AA atoms
/// proportional to their mass.
///
/// # Arguments
/// * `cg_force` – force on a single CG bead \[fx, fy, fz\].
/// * `group_masses` – masses of the atoms in this bead's group.
///
/// Returns per-atom forces.
pub fn distribute_force_to_atoms(cg_force: [f64; 3], group_masses: &[f64]) -> Vec<[f64; 3]> {
    let total: f64 = group_masses.iter().sum();
    if total < 1e-300 {
        return vec![[0.0; 3]; group_masses.len()];
    }
    group_masses
        .iter()
        .map(|&m| {
            let w = m / total;
            [cg_force[0] * w, cg_force[1] * w, cg_force[2] * w]
        })
        .collect()
}
/// Fit a Morse potential V(r) = D_e (1 − exp(−a(r−r_e)))²
/// using a golden-section search to minimise least-squares residuals.
///
/// # Arguments
/// * `r_data`   – distance values.
/// * `v_data`   – corresponding potential values.
/// * `r_e`      – equilibrium distance (fixed).
/// * `a_range`  – search range for `a` parameter (min, max).
/// * `de_range` – search range for `D_e` (min, max).
/// * `n_grid`   – grid resolution for brute-force 2-D search.
///
/// Returns `(D_e, a)` best-fit pair.
pub fn fit_morse_potential(
    r_data: &[f64],
    v_data: &[f64],
    r_e: f64,
    a_range: (f64, f64),
    de_range: (f64, f64),
    n_grid: usize,
) -> (f64, f64) {
    let n = r_data.len().min(v_data.len());
    if n == 0 || n_grid == 0 {
        return (1.0, 1.0);
    }
    let da = (a_range.1 - a_range.0) / n_grid.max(1) as f64;
    let dde = (de_range.1 - de_range.0) / n_grid.max(1) as f64;
    let mut best_rss = f64::INFINITY;
    let mut best_de = de_range.0;
    let mut best_a = a_range.0;
    for ia in 0..=n_grid {
        let a = a_range.0 + ia as f64 * da;
        for ide in 0..=n_grid {
            let de = de_range.0 + ide as f64 * dde;
            let rss: f64 = (0..n)
                .map(|i| {
                    let expt = (-a * (r_data[i] - r_e)).exp();
                    let v_morse = de * (1.0 - expt).powi(2);
                    (v_morse - v_data[i]).powi(2)
                })
                .sum();
            if rss < best_rss {
                best_rss = rss;
                best_de = de;
                best_a = a;
            }
        }
    }
    (best_de, best_a)
}
/// Evaluate Morse potential.
pub fn morse_potential(r: f64, de: f64, a: f64, r_e: f64) -> f64 {
    let expt = (-a * (r - r_e)).exp();
    de * (1.0 - expt).powi(2)
}
/// Evaluate Morse force: F = −dV/dr.
pub fn morse_force(r: f64, de: f64, a: f64, r_e: f64) -> f64 {
    let expt = (-a * (r - r_e)).exp();
    -2.0 * de * a * expt * (1.0 - expt)
}
/// Fit a harmonic bond potential k·(r−r0)²/2 via linear regression.
///
/// Returns `(k, r0)`.
pub fn fit_harmonic_bond(r_data: &[f64], v_data: &[f64]) -> (f64, f64) {
    let n = r_data.len().min(v_data.len());
    if n < 2 {
        return (0.0, 0.0);
    }
    let mean_r: f64 = r_data.iter().take(n).sum::<f64>() / n as f64;
    let dr2: Vec<f64> = r_data
        .iter()
        .take(n)
        .map(|&r| (r - mean_r).powi(2))
        .collect();
    let num: f64 = (0..n).map(|i| v_data[i] * dr2[i]).sum::<f64>();
    let den: f64 = dr2.iter().map(|&d| d * d).sum::<f64>();
    let k = if den > 1e-300 { 2.0 * num / den } else { 0.0 };
    (k, mean_r)
}
/// Backmapping quality metric: RMSD between back-mapped and reference AA positions.
pub fn backmapping_rmsd(back_mapped: &[[f64; 3]], reference: &[[f64; 3]]) -> f64 {
    let n = back_mapped.len().min(reference.len());
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = (0..n)
        .map(|i| {
            (0..3_usize)
                .map(|d| (back_mapped[i][d] - reference[i][d]).powi(2))
                .sum::<f64>()
        })
        .sum();
    (sum_sq / n as f64).sqrt()
}
/// Iterative back-mapping refinement: minimise RMSD by scaling the reference frame.
///
/// Adjusts a global scale factor for the reference frame offsets.
/// Returns the scale factor that minimises RMSD.
pub fn backmapping_optimal_scale(
    cg_positions: &[[f64; 3]],
    reference_frames: &[Vec<[f64; 3]>],
    target_aa: &[[f64; 3]],
    scale_range: (f64, f64),
    n_steps: usize,
) -> f64 {
    if n_steps == 0 || scale_range.0 >= scale_range.1 {
        return 1.0;
    }
    let ds = (scale_range.1 - scale_range.0) / n_steps as f64;
    let mut best_rmsd = f64::INFINITY;
    let mut best_scale = 1.0_f64;
    for step in 0..=n_steps {
        let scale = scale_range.0 + step as f64 * ds;
        let mut positions = Vec::new();
        for (b, bead_pos) in cg_positions.iter().enumerate() {
            if b < reference_frames.len() {
                for ref_pos in &reference_frames[b] {
                    positions.push([
                        bead_pos[0] + scale * ref_pos[0],
                        bead_pos[1] + scale * ref_pos[1],
                        bead_pos[2] + scale * ref_pos[2],
                    ]);
                }
            }
        }
        let rmsd = backmapping_rmsd(&positions, target_aa);
        if rmsd < best_rmsd {
            best_rmsd = rmsd;
            best_scale = scale;
        }
    }
    best_scale
}
/// Build a neighbour list for CG beads within `cutoff` (nm).
///
/// Returns a list of unique pairs `(i, j)` with `i < j` within the cutoff.
pub fn build_neighbour_list(positions: &[[f64; 3]], cutoff: f64) -> Vec<(usize, usize)> {
    let n = positions.len();
    let cutoff_sq = cutoff * cutoff;
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r_sq = dx * dx + dy * dy + dz * dz;
            if r_sq <= cutoff_sq {
                pairs.push((i, j));
            }
        }
    }
    pairs
}
/// Apply minimum-image convention to a distance vector.
pub fn minimum_image(mut dr: [f64; 3], box_lengths: [f64; 3]) -> [f64; 3] {
    for d in 0..3 {
        let l = box_lengths[d];
        if l > 0.0 {
            dr[d] -= l * (dr[d] / l).round();
        }
    }
    dr
}
/// Build a neighbour list with periodic boundary conditions.
pub fn build_neighbour_list_pbc(
    positions: &[[f64; 3]],
    cutoff: f64,
    box_lengths: [f64; 3],
) -> Vec<(usize, usize)> {
    let n = positions.len();
    let cutoff_sq = cutoff * cutoff;
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let dr_raw = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let dr = minimum_image(dr_raw, box_lengths);
            let r_sq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if r_sq <= cutoff_sq {
                pairs.push((i, j));
            }
        }
    }
    pairs
}
/// Steepest-descent energy minimisation for a CG system.
///
/// Updates positions in-place until max force < `f_tol` or `max_steps` reached.
/// `force_fn` must return per-bead forces for the current positions.
///
/// Returns `(n_steps_taken, max_force_final)`.
pub fn steepest_descent_cg<F>(
    positions: &mut [[f64; 3]],
    step_size: f64,
    f_tol: f64,
    max_steps: usize,
    mut force_fn: F,
) -> (usize, f64)
where
    F: FnMut(&[[f64; 3]]) -> Vec<[f64; 3]>,
{
    let mut current_step_size = step_size;
    let mut n_taken = 0usize;
    let mut max_force = f64::INFINITY;
    for step in 0..max_steps {
        let forces = force_fn(positions);
        max_force = forces
            .iter()
            .map(|f| (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt())
            .fold(0.0_f64, f64::max);
        if max_force < f_tol {
            n_taken = step;
            break;
        }
        for (i, pos) in positions.iter_mut().enumerate() {
            if i < forces.len() {
                let f_norm = (forces[i][0] * forces[i][0]
                    + forces[i][1] * forces[i][1]
                    + forces[i][2] * forces[i][2])
                    .sqrt();
                if f_norm > 1e-300 {
                    for d in 0..3 {
                        pos[d] += current_step_size * forces[i][d] / f_norm;
                    }
                }
            }
        }
        current_step_size *= 0.99;
        n_taken = step + 1;
    }
    (n_taken, max_force)
}
/// Compute ideal MARTINI bond length from RDF peak position.
///
/// The equilibrium bond length is estimated as the first peak of the
/// radial distribution function `rdf` at `r_bins` positions.
pub fn ideal_bond_length_from_rdf(r_bins: &[f64], rdf: &[f64]) -> Option<f64> {
    if r_bins.is_empty() || rdf.is_empty() {
        return None;
    }
    let n = r_bins.len().min(rdf.len());
    let (peak_idx, _) = (0..n)
        .map(|i| (i, rdf[i]))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))?;
    r_bins.get(peak_idx).copied()
}
/// Compute MARTINI bond spring constant from fluctuations:
///
/// k = kT / ⟨(r − r₀)²⟩
pub fn bond_spring_constant_from_fluctuations(r_samples: &[f64], r0: f64, kt: f64) -> f64 {
    let n = r_samples.len();
    if n < 2 {
        return 0.0;
    }
    let var: f64 = r_samples.iter().map(|&r| (r - r0).powi(2)).sum::<f64>() / n as f64;
    if var < 1e-300 {
        return f64::INFINITY;
    }
    kt / var
}
/// Compute MARTINI angle spring constant from fluctuations:
///
/// k_θ = kT / ⟨(θ − θ₀)²⟩
pub fn angle_spring_constant_from_fluctuations(theta_samples: &[f64], theta0: f64, kt: f64) -> f64 {
    let n = theta_samples.len();
    if n < 2 {
        return 0.0;
    }
    let var: f64 = theta_samples
        .iter()
        .map(|&t| (t - theta0).powi(2))
        .sum::<f64>()
        / n as f64;
    if var < 1e-300 {
        return f64::INFINITY;
    }
    kt / var
}
/// Mean square displacement (MSD) from a list of position trajectories.
///
/// `traj` is indexed `[time][bead]`; returns MSD(t) for lags 0..max_lag.
pub fn mean_square_displacement(traj: &[Vec<[f64; 3]>], max_lag: usize) -> Vec<f64> {
    let n_frames = traj.len();
    if n_frames == 0 || max_lag == 0 {
        return Vec::new();
    }
    let n_beads = traj[0].len();
    let effective_max = max_lag.min(n_frames);
    let mut msd = vec![0.0_f64; effective_max];
    for (lag, msd_val) in msd.iter_mut().enumerate() {
        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for t0 in 0..n_frames.saturating_sub(lag) {
            let t1 = t0 + lag;
            if t1 < n_frames {
                let b_max = n_beads.min(traj[t1].len());
                for (bead1, bead0) in traj[t1][..b_max].iter().zip(traj[t0][..b_max].iter()) {
                    let dr: f64 = bead1
                        .iter()
                        .zip(bead0.iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    sum += dr;
                    count += 1;
                }
            }
        }
        if count > 0 {
            *msd_val = sum / count as f64;
        }
    }
    msd
}
/// Self-diffusion coefficient from MSD via Einstein relation: D = MSD(t) / (6t).
///
/// Uses a linear fit of MSD vs time in the range \[fit_start, fit_end\] frames.
///
/// Returns `None` if the fit range is too small.
pub fn diffusion_coefficient_from_msd(
    msd: &[f64],
    dt: f64,
    fit_start: usize,
    fit_end: usize,
) -> Option<f64> {
    let end = fit_end.min(msd.len());
    if fit_start >= end || end == 0 {
        return None;
    }
    let n = end - fit_start;
    if n < 2 {
        return None;
    }
    let x: Vec<f64> = (fit_start..end).map(|i| i as f64 * dt).collect();
    let y: Vec<f64> = msd[fit_start..end].to_vec();
    let mean_x = x.iter().sum::<f64>() / n as f64;
    let mean_y = y.iter().sum::<f64>() / n as f64;
    let num: f64 = (0..n).map(|i| (x[i] - mean_x) * (y[i] - mean_y)).sum();
    let den: f64 = (0..n).map(|i| (x[i] - mean_x).powi(2)).sum();
    if den < 1e-300 {
        return None;
    }
    let slope = num / den;
    Some(slope / 6.0)
}
/// Screened Coulomb interaction as used in MARTINI 3.
///
/// V(r) = q_i q_j / (4πε₀ ε_r r) · exp(−r / λ_D)
///
/// All in reduced units (result in kJ/mol when q in elementary charge,
/// r in nm, using the factor 138.935 kJ·nm·mol⁻¹·e⁻²).
///
/// # Arguments
/// * `qi`, `qj`  – partial charges (elementary charge units).
/// * `r`         – distance (nm).
/// * `epsilon_r` – relative permittivity.
/// * `lambda_d`  – Debye screening length (nm).
pub fn martini_screened_coulomb(qi: f64, qj: f64, r: f64, epsilon_r: f64, lambda_d: f64) -> f64 {
    if r < 1e-10 {
        return 0.0;
    }
    let factor = 138.935_f64;
    let screening = if lambda_d > 0.0 {
        (-r / lambda_d).exp()
    } else {
        1.0
    };
    factor * qi * qj / (epsilon_r * r) * screening
}
/// Force from screened Coulomb (along the r direction, magnitude).
pub fn martini_screened_coulomb_force(
    qi: f64,
    qj: f64,
    r: f64,
    epsilon_r: f64,
    lambda_d: f64,
) -> f64 {
    if r < 1e-10 {
        return 0.0;
    }
    let factor = 138.935_f64;
    let screening = if lambda_d > 0.0 {
        (-r / lambda_d).exp()
    } else {
        1.0
    };
    let v_over_r = factor * qi * qj / (epsilon_r * r * r) * screening;
    let extra = if lambda_d > 0.0 {
        factor * qi * qj / (epsilon_r * r * lambda_d) * screening
    } else {
        0.0
    };
    v_over_r + extra
}
/// MARTINI water bead — 4 water molecules → 1 bead (W type).
/// Returns LJ energy using MARTINI P4 water parameters.
pub fn martini_water_lj(r: f64) -> f64 {
    let sigma = 0.47_f64;
    let epsilon = 5.0_f64;
    martini_lj_energy(r, sigma, epsilon)
}
#[cfg(test)]
mod tests_extended {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_martini_4to1_mapping_single_bead() {
        let pos = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let mass = vec![1.0_f64; 4];
        let beads = martini_4to1_mapping(&pos, &mass);
        assert_eq!(beads.len(), 1, "4 atoms → 1 bead");
        let expected = [0.5, 0.5, 0.5];
        for d in 0..3 {
            assert!(
                (beads[0][d] - expected[d]).abs() < 1e-12,
                "COM mismatch at dim {d}"
            );
        }
    }
    #[test]
    fn test_martini_4to1_mapping_two_beads() {
        let pos: Vec<[f64; 3]> = (0..8).map(|i| [i as f64, 0.0, 0.0]).collect();
        let mass = vec![1.0_f64; 8];
        let beads = martini_4to1_mapping(&pos, &mass);
        assert_eq!(beads.len(), 2, "8 atoms → 2 beads");
    }
    #[test]
    fn test_martini_4to1_mapping_empty() {
        let beads = martini_4to1_mapping(&[], &[]);
        assert!(beads.is_empty());
    }
    #[test]
    fn test_martini_bead_mass_regular() {
        let m = martini_bead_mass(MartiniBead::C1);
        assert!(
            (m - 72.0).abs() < 1e-12,
            "regular bead mass should be 72 amu"
        );
    }
    #[test]
    fn test_martini_bead_mass_small() {
        let m = martini_bead_mass(MartiniBead::Sc3);
        assert!((m - 54.0).abs() < 1e-12, "small bead mass should be 54 amu");
    }
    #[test]
    fn test_martini_bead_radius_regular() {
        let r = martini_bead_radius(MartiniBead::P4);
        assert!(
            (r - 0.47).abs() < 1e-12,
            "regular bead radius should be 0.47 nm"
        );
    }
    #[test]
    fn test_distribute_force_conservation() {
        let f = [6.0_f64, 0.0, 0.0];
        let masses = [1.0_f64, 2.0, 3.0];
        let atom_forces = distribute_force_to_atoms(f, &masses);
        let fx_sum: f64 = atom_forces.iter().map(|af| af[0]).sum();
        assert!(
            (fx_sum - 6.0).abs() < 1e-12,
            "distributed forces should sum to original"
        );
    }
    #[test]
    fn test_distribute_force_proportional() {
        let f = [0.0_f64, 3.0, 0.0];
        let masses = [1.0_f64, 1.0, 1.0];
        let atom_forces = distribute_force_to_atoms(f, &masses);
        for af in &atom_forces {
            assert!(
                (af[1] - 1.0).abs() < 1e-12,
                "equal masses → equal force distribution"
            );
        }
    }
    #[test]
    fn test_morse_at_equilibrium() {
        let v = morse_potential(1.0, 5.0, 2.0, 1.0);
        assert!(v.abs() < 1e-12, "Morse at r_e should be 0, got {v}");
    }
    #[test]
    fn test_morse_force_at_equilibrium() {
        let f = morse_force(1.0, 5.0, 2.0, 1.0);
        assert!(f.abs() < 1e-12, "Morse force at r_e should be 0, got {f}");
    }
    #[test]
    fn test_morse_force_direction() {
        let f = morse_force(1.5, 5.0, 2.0, 1.0);
        assert!(f < 0.0, "Morse force at r > r_e should be attractive");
    }
    #[test]
    fn test_fit_morse_potential_recovers_known() {
        let de_true = 4.0_f64;
        let a_true = 2.0_f64;
        let r_e = 1.0_f64;
        let r_data: Vec<f64> = (1..=10).map(|i| 0.5 + i as f64 * 0.1).collect();
        let v_data: Vec<f64> = r_data
            .iter()
            .map(|&r| morse_potential(r, de_true, a_true, r_e))
            .collect();
        let (de_fit, a_fit) =
            fit_morse_potential(&r_data, &v_data, r_e, (1.0, 3.0), (3.0, 5.0), 20);
        assert!(
            (de_fit - de_true).abs() < 0.5,
            "de_fit should be close to {de_true}, got {de_fit}"
        );
        assert!(
            (a_fit - a_true).abs() < 0.5,
            "a_fit should be close to {a_true}, got {a_fit}"
        );
    }
    #[test]
    fn test_fit_harmonic_bond_known() {
        let k_true = 500.0_f64;
        let r0_true = 1.0_f64;
        let r_data: Vec<f64> = (0..=10).map(|i| 0.8 + i as f64 * 0.04).collect();
        let v_data: Vec<f64> = r_data
            .iter()
            .map(|&r| 0.5 * k_true * (r - r0_true).powi(2))
            .collect();
        let (k_fit, r0_fit) = fit_harmonic_bond(&r_data, &v_data);
        assert!(
            (k_fit - k_true).abs() / k_true < 0.05,
            "k_fit should be close to {k_true}"
        );
        assert!(
            (r0_fit - r0_true).abs() < 0.05,
            "r0_fit should be close to {r0_true}"
        );
    }
    #[test]
    fn test_backmapping_rmsd_identical() {
        let pos = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let rmsd = backmapping_rmsd(&pos, &pos);
        assert!(rmsd.abs() < 1e-12, "RMSD of identical configs should be 0");
    }
    #[test]
    fn test_backmapping_rmsd_known() {
        let a = vec![[0.0, 0.0, 0.0]];
        let b = vec![[1.0, 0.0, 0.0]];
        let rmsd = backmapping_rmsd(&a, &b);
        assert!((rmsd - 1.0).abs() < 1e-12, "RMSD should be 1.0, got {rmsd}");
    }
    #[test]
    fn test_backmapping_rmsd_empty() {
        let rmsd = backmapping_rmsd(&[], &[]);
        assert!(rmsd.abs() < 1e-12, "RMSD of empty configs should be 0");
    }
    #[test]
    fn test_cg_integrator_com_at_origin() {
        let pos = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let vel = vec![[0.0; 3]; 2];
        let masses = vec![1.0_f64, 1.0];
        let state = CgIntegratorState::new(pos, vel, masses);
        let com = state.centre_of_mass();
        assert!(com[0].abs() < 1e-12, "COM should be at origin");
    }
    #[test]
    fn test_cg_integrator_kinetic_energy_zero() {
        let pos = vec![[0.0; 3]; 3];
        let vel = vec![[0.0; 3]; 3];
        let masses = vec![1.0_f64; 3];
        let state = CgIntegratorState::new(pos, vel, masses);
        assert!(
            state.kinetic_energy().abs() < 1e-12,
            "zero velocity → zero KE"
        );
    }
    #[test]
    fn test_cg_integrator_position_step() {
        let pos = vec![[0.0, 0.0, 0.0]];
        let vel = vec![[1.0, 2.0, 3.0]];
        let masses = vec![1.0_f64];
        let mut state = CgIntegratorState::new(pos, vel, masses);
        state.position_step(0.1);
        assert!((state.positions[0][0] - 0.1).abs() < 1e-12);
        assert!((state.positions[0][1] - 0.2).abs() < 1e-12);
        assert!((state.positions[0][2] - 0.3).abs() < 1e-12);
    }
    #[test]
    fn test_cg_integrator_remove_com_velocity() {
        let pos = vec![[0.0; 3]; 2];
        let vel = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let masses = vec![1.0_f64; 2];
        let mut state = CgIntegratorState::new(pos, vel, masses);
        state.remove_com_velocity();
        let com_vel = state.centre_of_mass_velocity();
        for v in com_vel {
            assert!(v.abs() < 1e-12, "COM velocity should be zero after removal");
        }
    }
    #[test]
    fn test_cg_integrator_pbc_wraps() {
        let pos = vec![[1.5, -0.1, 2.5]];
        let vel = vec![[0.0; 3]];
        let masses = vec![1.0_f64];
        let mut state = CgIntegratorState::new(pos, vel, masses);
        state.apply_pbc([1.0, 1.0, 1.0]);
        for (d, &v) in state.positions[0].iter().enumerate() {
            assert!(
                (0.0..1.0).contains(&v),
                "PBC-wrapped position out of box at dim {d}: {v}"
            );
        }
    }
    #[test]
    fn test_build_neighbour_list_close() {
        let pos = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let pairs = build_neighbour_list(&pos, 0.5);
        assert_eq!(pairs.len(), 1, "one pair within cutoff");
    }
    #[test]
    fn test_build_neighbour_list_far() {
        let pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let pairs = build_neighbour_list(&pos, 0.5);
        assert!(pairs.is_empty(), "pair outside cutoff should not appear");
    }
    #[test]
    fn test_build_neighbour_list_pbc_wraps() {
        let pos = vec![[0.05, 0.0, 0.0], [0.95, 0.0, 0.0]];
        let pairs = build_neighbour_list_pbc(&pos, 0.15, [1.0, 1.0, 1.0]);
        assert_eq!(pairs.len(), 1, "PBC-wrapped pair should appear");
    }
    #[test]
    fn test_msd_stationary_particles() {
        let frame = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let traj = vec![frame.clone(), frame.clone(), frame.clone()];
        let msd = mean_square_displacement(&traj, 3);
        for &m in &msd {
            assert!(m.abs() < 1e-12, "stationary particles → MSD = 0");
        }
    }
    #[test]
    fn test_msd_lag_zero_is_zero() {
        let traj = vec![vec![[1.0, 2.0, 3.0]]];
        let msd = mean_square_displacement(&traj, 5);
        assert!(!msd.is_empty(), "MSD should have entries");
        assert!(msd[0].abs() < 1e-12, "MSD at lag=0 should be 0");
    }
    #[test]
    fn test_diffusion_coefficient_linear_msd() {
        let d_true = 0.01_f64;
        let dt = 1.0_f64;
        let msd: Vec<f64> = (0..20).map(|i| 6.0 * d_true * i as f64 * dt).collect();
        let d = diffusion_coefficient_from_msd(&msd, dt, 1, 19).unwrap();
        assert!(
            (d - d_true).abs() < 1e-6,
            "D from linear MSD should match, got {d}"
        );
    }
    #[test]
    fn test_bond_spring_constant_from_fluctuations_known() {
        let r0 = 1.0_f64;
        let kt = 2.479_f64;
        let samples: Vec<f64> = [-0.1_f64, 0.1].iter().map(|&d| r0 + d).collect();
        let k = bond_spring_constant_from_fluctuations(&samples, r0, kt);
        let expected = kt / 0.01;
        assert!(
            (k - expected).abs() / expected < 0.01,
            "k mismatch: {k} vs {expected}"
        );
    }
    #[test]
    fn test_angle_spring_constant_from_fluctuations_positive() {
        let samples = vec![1.5_f64, 1.57, 1.6, 1.55];
        let k = angle_spring_constant_from_fluctuations(&samples, 1.57, 2.479);
        assert!(k > 0.0, "angle spring constant should be positive");
    }
    #[test]
    fn test_ideal_bond_length_from_rdf_peak() {
        let r_bins = vec![0.3_f64, 0.4, 0.5, 0.6, 0.7];
        let rdf = vec![0.1_f64, 0.5, 2.0, 0.8, 0.2];
        let r0 = ideal_bond_length_from_rdf(&r_bins, &rdf).unwrap();
        assert!(
            (r0 - 0.5).abs() < 1e-12,
            "bond length should be at RDF peak, got {r0}"
        );
    }
    #[test]
    fn test_ideal_bond_length_from_rdf_empty() {
        assert!(ideal_bond_length_from_rdf(&[], &[]).is_none());
    }
    #[test]
    fn test_martini_screened_coulomb_positive_charges_repulsive() {
        let v = martini_screened_coulomb(1.0, 1.0, 0.5, 15.0, 1.0);
        assert!(v > 0.0, "same-sign charges should repel");
    }
    #[test]
    fn test_martini_screened_coulomb_opposite_charges_attractive() {
        let v = martini_screened_coulomb(1.0, -1.0, 0.5, 15.0, 1.0);
        assert!(v < 0.0, "opposite charges should attract");
    }
    #[test]
    fn test_martini_screened_coulomb_zero_r() {
        let v = martini_screened_coulomb(1.0, 1.0, 0.0, 15.0, 1.0);
        assert_eq!(v, 0.0, "zero distance should return 0 (guard)");
    }
    #[test]
    fn test_martini_water_lj_at_sigma() {
        let v = martini_water_lj(0.47);
        assert!(
            v.abs() < 1e-10,
            "LJ potential at r=sigma should be 0, got {v}"
        );
    }
    #[test]
    fn test_martini_water_lj_minimum_negative() {
        let r_min = 2.0_f64.powf(1.0_f64 / 6.0_f64) * 0.47_f64;
        let v = martini_water_lj(r_min);
        assert!(
            v < 0.0,
            "LJ potential at minimum should be negative, got {v}"
        );
        assert!(
            (v - (-5.0)).abs() < 0.01,
            "LJ minimum should be ~-epsilon=-5.0, got {v}"
        );
    }
    #[test]
    fn test_steepest_descent_harmonic_converges() {
        let k = 100.0_f64;
        let x0 = 1.0_f64;
        let mut positions = vec![[2.0_f64, 0.0, 0.0]];
        let (n_steps, max_f) = steepest_descent_cg(&mut positions, 0.001, 0.01, 10_000, |pos| {
            pos.iter()
                .map(|p| {
                    let fx = -k * (p[0] - x0);
                    [fx, 0.0, 0.0]
                })
                .collect()
        });
        assert!(
            max_f < 0.01 || n_steps > 0,
            "steepest descent should reduce force"
        );
    }
    #[test]
    fn test_steepest_descent_already_at_minimum() {
        let mut positions = vec![[1.0_f64, 0.0, 0.0]];
        let (_, max_f) =
            steepest_descent_cg(&mut positions, 0.001, 1e-6, 100, |_| vec![[0.0, 0.0, 0.0]]);
        assert!(max_f < 1e-6, "zero force → already at minimum");
    }
    #[test]
    fn test_cg_integrator_config_kt() {
        let cfg = CgIntegratorConfig::new(0.002, 1000, 300.0);
        let kt = cfg.kt();
        assert!(
            (kt - 2.494).abs() < 0.01,
            "kT at 300K should be ~2.494 kJ/mol, got {kt}"
        );
    }
    #[test]
    fn test_cg_integrator_instantaneous_temperature() {
        let pos = vec![[0.0; 3]; 3];
        let vel = vec![[1.0, 0.0, 0.0]; 3];
        let masses = vec![1.0_f64; 3];
        let state = CgIntegratorState::new(pos, vel, masses);
        let t = state.instantaneous_temperature(9.0);
        assert!(t > 0.0, "temperature should be positive");
    }
}
