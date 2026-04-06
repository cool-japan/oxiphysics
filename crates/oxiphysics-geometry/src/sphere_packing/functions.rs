//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PackingSphere;

/// Compute the packing fraction of a set of spheres in a given domain volume.
pub fn compute_packing_fraction(spheres: &[PackingSphere], domain_volume: f64) -> f64 {
    let sv: f64 = spheres.iter().map(|s| s.volume()).sum();
    if domain_volume > 1e-30 {
        sv / domain_volume
    } else {
        0.0
    }
}
/// Compute the coordination number distribution.
///
/// Returns a vector of (coordination_number, count) pairs.
pub fn coordination_number_distribution(
    spheres: &[PackingSphere],
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let n = spheres.len();
    let mut coord_nums = vec![0usize; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let ov = spheres[i].overlap_with(&spheres[j]);
            if ov.abs() < tolerance {
                coord_nums[i] += 1;
                coord_nums[j] += 1;
            }
        }
    }
    let max_cn = *coord_nums.iter().max().unwrap_or(&0);
    let mut dist = vec![(0usize, 0usize); max_cn + 1];
    for (i, cn) in dist.iter_mut().enumerate() {
        *cn = (i, 0);
    }
    for &cn in &coord_nums {
        dist[cn].1 += 1;
    }
    dist.retain(|(_, count)| *count > 0);
    dist
}
/// Check if any two spheres in the list overlap.
pub fn has_overlap(spheres: &[PackingSphere]) -> bool {
    for i in 0..spheres.len() {
        for j in (i + 1)..spheres.len() {
            if spheres[i].overlaps(&spheres[j]) {
                return true;
            }
        }
    }
    false
}
/// Compute the second-order structure factor S(q) for a set of sphere centers.
///
/// `q_values` are the wavenumbers to evaluate at.
pub fn structure_factor(spheres: &[PackingSphere], q_values: &[f64]) -> Vec<f64> {
    let n = spheres.len() as f64;
    if n < 1.0 {
        return vec![0.0; q_values.len()];
    }
    q_values
        .iter()
        .map(|&q| {
            let mut sum_cos = 0.0f64;
            let mut sum_sin = 0.0f64;
            for s in spheres {
                let phase = q * (s.center[0] + s.center[1] + s.center[2]);
                sum_cos += phase.cos();
                sum_sin += phase.sin();
            }
            (sum_cos * sum_cos + sum_sin * sum_sin) / n
        })
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sphere_packing::ConfinementConfig;
    use crate::sphere_packing::ConfinementShape;
    use crate::sphere_packing::ConfiningGeometry;
    use crate::sphere_packing::LatticeType;
    use crate::sphere_packing::OptimizerConfig;
    use crate::sphere_packing::OrderedPacking;
    use crate::sphere_packing::PackingOptimizer;
    use crate::sphere_packing::PolyDispersePacking;
    use crate::sphere_packing::RandomPacking;
    use crate::sphere_packing::RsaConfig;
    use crate::sphere_packing::SizeDistribution;
    use crate::sphere_packing::SizeDistributionParams;
    use crate::sphere_packing::VoidSpaceAnalysis;
    #[test]
    fn test_sphere_volume() {
        let s = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let expected = 4.0 / 3.0 * std::f64::consts::PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_sphere_overlap_true() {
        let s1 = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let s2 = PackingSphere::new([1.5, 0.0, 0.0], 1.0, 1);
        assert!(s1.overlaps(&s2));
    }
    #[test]
    fn test_sphere_overlap_false() {
        let s1 = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let s2 = PackingSphere::new([3.0, 0.0, 0.0], 1.0, 1);
        assert!(!s1.overlaps(&s2));
    }
    #[test]
    fn test_sphere_distance() {
        let s1 = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let s2 = PackingSphere::new([3.0, 4.0, 0.0], 1.0, 1);
        assert!((s1.distance_to(&s2) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_sphere_overlap_with() {
        let s1 = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let s2 = PackingSphere::new([1.5, 0.0, 0.0], 1.0, 1);
        let ov = s1.overlap_with(&s2);
        assert!((ov - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_rsa_creation() {
        let config = RsaConfig {
            radius: 0.1,
            domain: [1.0, 1.0, 1.0],
            max_attempts: 100,
            min_gap: 0.0,
        };
        let rsa = RandomPacking::new(config);
        assert_eq!(rsa.spheres.len(), 0);
        assert!(!rsa.jammed);
    }
    #[test]
    fn test_rsa_place_at_no_overlap() {
        let config = RsaConfig {
            radius: 0.1,
            domain: [2.0, 2.0, 2.0],
            max_attempts: 100,
            min_gap: 0.0,
        };
        let mut rsa = RandomPacking::new(config);
        assert!(rsa.place_at([0.5, 0.5, 0.5]));
        assert!(rsa.place_at([1.5, 0.5, 0.5]));
        assert_eq!(rsa.spheres.len(), 2);
    }
    #[test]
    fn test_rsa_place_at_overlap_rejected() {
        let config = RsaConfig {
            radius: 0.2,
            domain: [2.0, 2.0, 2.0],
            max_attempts: 100,
            min_gap: 0.0,
        };
        let mut rsa = RandomPacking::new(config);
        assert!(rsa.place_at([1.0, 1.0, 1.0]));
        assert!(!rsa.place_at([1.1, 1.0, 1.0]));
    }
    #[test]
    fn test_rsa_packing_fraction_positive() {
        let config = RsaConfig {
            radius: 0.05,
            domain: [1.0, 1.0, 1.0],
            max_attempts: 50,
            min_gap: 0.0,
        };
        let mut rsa = RandomPacking::new(config);
        rsa.run();
        let pf = rsa.packing_fraction();
        assert!(pf > 0.0 && pf <= 1.0);
    }
    #[test]
    fn test_rsa_jamming_flag() {
        let config = RsaConfig {
            radius: 0.4,
            domain: [1.0, 1.0, 1.0],
            max_attempts: 20,
            min_gap: 0.0,
        };
        let mut rsa = RandomPacking::new(config);
        rsa.run();
        assert!(rsa.jammed);
    }
    #[test]
    fn test_rsa_theoretical_jamming_limit() {
        let limit = RandomPacking::theoretical_jamming_limit();
        assert!((limit - 0.3841).abs() < 1e-4);
    }
    #[test]
    fn test_rsa_valid_placement() {
        let config = RsaConfig {
            radius: 0.1,
            domain: [1.0, 1.0, 1.0],
            max_attempts: 100,
            min_gap: 0.0,
        };
        let mut rsa = RandomPacking::new(config);
        rsa.place_at([0.5, 0.5, 0.5]);
        let c = PackingSphere::new([0.8, 0.5, 0.5], 0.1, 1);
        assert!(rsa.is_valid_placement(&c));
        let c2 = PackingSphere::new([0.55, 0.5, 0.5], 0.1, 2);
        assert!(!rsa.is_valid_placement(&c2));
    }
    #[test]
    fn test_ordered_fcc_generation() {
        let mut packing = OrderedPacking::new(LatticeType::Fcc, 1.0, [2, 2, 2]);
        packing.generate();
        assert_eq!(packing.spheres.len(), 32);
    }
    #[test]
    fn test_ordered_sc_generation() {
        let mut packing = OrderedPacking::new(LatticeType::SimpleCubic, 1.0, [3, 3, 3]);
        packing.generate();
        assert_eq!(packing.spheres.len(), 27);
    }
    #[test]
    fn test_ordered_bcc_generation() {
        let mut packing = OrderedPacking::new(LatticeType::Bcc, 1.0, [2, 2, 2]);
        packing.generate();
        assert_eq!(packing.spheres.len(), 16);
    }
    #[test]
    fn test_fcc_theoretical_packing_fraction() {
        let packing = OrderedPacking::new(LatticeType::Fcc, 1.0, [1, 1, 1]);
        let eta = packing.theoretical_packing_fraction();
        assert!((eta - 0.7405).abs() < 0.0001);
    }
    #[test]
    fn test_bcc_theoretical_packing_fraction() {
        let packing = OrderedPacking::new(LatticeType::Bcc, 1.0, [1, 1, 1]);
        let eta = packing.theoretical_packing_fraction();
        assert!((eta - 0.6802).abs() < 0.0001);
    }
    #[test]
    fn test_sc_theoretical_packing_fraction() {
        let packing = OrderedPacking::new(LatticeType::SimpleCubic, 1.0, [1, 1, 1]);
        let eta = packing.theoretical_packing_fraction();
        assert!((eta - std::f64::consts::FRAC_PI_6).abs() < 0.0001);
    }
    #[test]
    fn test_fcc_coordination_number() {
        let packing = OrderedPacking::new(LatticeType::Fcc, 1.0, [1, 1, 1]);
        assert_eq!(packing.coordination_number(), 12);
    }
    #[test]
    fn test_bcc_coordination_number() {
        let packing = OrderedPacking::new(LatticeType::Bcc, 1.0, [1, 1, 1]);
        assert_eq!(packing.coordination_number(), 8);
    }
    #[test]
    fn test_ordered_no_overlap() {
        let mut packing = OrderedPacking::new(LatticeType::Fcc, 1.0, [2, 2, 2]);
        packing.generate();
        assert!(!has_overlap(&packing.spheres));
    }
    #[test]
    fn test_ordered_find_neighbors() {
        let mut packing = OrderedPacking::new(LatticeType::SimpleCubic, 1.0, [3, 3, 3]);
        packing.generate();
        let a = packing.lattice_parameter;
        let neighbors = packing.find_neighbors(13, a * 1.01);
        assert!(neighbors.len() >= 4);
    }
    #[test]
    fn test_nearest_neighbor_distance() {
        let packing = OrderedPacking::new(LatticeType::Fcc, 2.5, [1, 1, 1]);
        assert!((packing.nearest_neighbor_distance() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_polydisperse_monodisperse() {
        let params = SizeDistributionParams {
            mean_radius: 0.1,
            ..Default::default()
        };
        let mut packing = PolyDispersePacking::new(
            SizeDistribution::Monodisperse,
            params,
            [2.0, 2.0, 2.0],
            1000,
        );
        packing.run(10);
        for s in &packing.spheres {
            assert!((s.radius - 0.1).abs() < 1e-12);
        }
    }
    #[test]
    fn test_polydisperse_bidisperse() {
        let params = SizeDistributionParams {
            mean_radius: 0.05,
            size_ratio: 2.0,
            large_fraction: 0.5,
            ..Default::default()
        };
        let mut packing =
            PolyDispersePacking::new(SizeDistribution::Bidisperse, params, [2.0, 2.0, 2.0], 2000);
        packing.run(20);
        let radii: std::collections::HashSet<_> = packing
            .spheres
            .iter()
            .map(|s| (s.radius * 1e10).round() as i64)
            .collect();
        assert!(radii.len() <= 2);
    }
    #[test]
    fn test_polydisperse_no_overlap() {
        let params = SizeDistributionParams::default();
        let mut packing = PolyDispersePacking::new(
            SizeDistribution::Monodisperse,
            params,
            [3.0, 3.0, 3.0],
            1000,
        );
        packing.run(10);
        assert!(!has_overlap(&packing.spheres));
    }
    #[test]
    fn test_polydisperse_packing_fraction_positive() {
        let params = SizeDistributionParams::default();
        let mut packing = PolyDispersePacking::new(
            SizeDistribution::Monodisperse,
            params,
            [2.0, 2.0, 2.0],
            1000,
        );
        packing.run(5);
        assert!(packing.packing_fraction() > 0.0);
    }
    #[test]
    fn test_polydisperse_radius_statistics() {
        let params = SizeDistributionParams {
            mean_radius: 0.05,
            ..Default::default()
        };
        let mut packing = PolyDispersePacking::new(
            SizeDistribution::Monodisperse,
            params,
            [2.0, 2.0, 2.0],
            1000,
        );
        packing.run(5);
        let (mean, std, min, max) = packing.radius_statistics();
        assert!(mean > 0.0);
        assert!(min <= mean);
        assert!(max >= mean);
        let _ = std;
    }
    #[test]
    fn test_polydisperse_pdi_monodisperse() {
        let params = SizeDistributionParams {
            mean_radius: 0.05,
            ..Default::default()
        };
        let mut packing = PolyDispersePacking::new(
            SizeDistribution::Monodisperse,
            params,
            [2.0, 2.0, 2.0],
            1000,
        );
        packing.run(5);
        let pdi = packing.polydispersity_index();
        assert!(pdi < 1e-10, "PDI for monodisperse should be 0, got {}", pdi);
    }
    #[test]
    fn test_polydisperse_histogram() {
        let params = SizeDistributionParams {
            mean_radius: 0.05,
            ..Default::default()
        };
        let mut packing = PolyDispersePacking::new(
            SizeDistribution::Monodisperse,
            params,
            [2.0, 2.0, 2.0],
            1000,
        );
        packing.run(5);
        let (bins, counts) = packing.size_histogram(5);
        let total: usize = counts.iter().sum();
        assert_eq!(total, packing.spheres.len());
        assert_eq!(bins.len(), counts.len());
    }
    #[test]
    fn test_sample_radius_lognormal_bounds() {
        let params = SizeDistributionParams {
            mean_radius: 0.05,
            std_radius: 0.2,
            min_radius: 0.01,
            max_radius: 0.1,
            ..Default::default()
        };
        let packing = PolyDispersePacking::new(
            SizeDistribution::LogNormal,
            params.clone(),
            [2.0, 2.0, 2.0],
            1000,
        );
        let mut rng = rand::rng();
        for _ in 0..100 {
            let r = packing.sample_radius(&mut rng);
            assert!(r >= params.min_radius && r <= params.max_radius);
        }
    }
    #[test]
    fn test_optimizer_creation() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([1.5, 0.0, 0.0], 1.0, 1),
        ];
        let opt = PackingOptimizer::new(spheres, [5.0, 5.0, 5.0], OptimizerConfig::default());
        assert_eq!(opt.spheres.len(), 2);
    }
    #[test]
    fn test_optimizer_no_overlap_zero_energy() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([3.0, 0.0, 0.0], 1.0, 1),
        ];
        let opt = PackingOptimizer::new(spheres, [10.0, 10.0, 10.0], OptimizerConfig::default());
        assert_eq!(opt.potential_energy(), 0.0);
    }
    #[test]
    fn test_optimizer_overlap_nonzero_energy() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([1.5, 0.0, 0.0], 1.0, 1),
        ];
        let opt = PackingOptimizer::new(spheres, [10.0, 10.0, 10.0], OptimizerConfig::default());
        assert!(opt.potential_energy() > 0.0);
    }
    #[test]
    fn test_optimizer_step_reduces_overlap() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([1.5, 0.0, 0.0], 1.0, 1),
        ];
        let mut opt =
            PackingOptimizer::new(spheres, [10.0, 10.0, 10.0], OptimizerConfig::default());
        let before = opt.max_overlap();
        opt.step();
        let after = opt.max_overlap();
        assert!(
            after <= before,
            "overlap should decrease: before={}, after={}",
            before,
            after
        );
    }
    #[test]
    fn test_optimizer_run_converges_no_overlap() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([1.5, 0.0, 0.0], 1.0, 1),
        ];
        let config = OptimizerConfig {
            max_iterations: 5000,
            tolerance: 1e-5,
            time_step: 0.05,
            damping: 0.3,
            stiffness: 1.0,
        };
        let mut opt = PackingOptimizer::new(spheres, [10.0, 10.0, 10.0], config);
        opt.run();
        assert!(opt.converged || opt.max_overlap() < 1e-4);
    }
    #[test]
    fn test_optimizer_packing_fraction() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 0.5, 0),
            PackingSphere::new([2.0, 0.0, 0.0], 0.5, 1),
        ];
        let opt = PackingOptimizer::new(spheres, [4.0, 4.0, 4.0], OptimizerConfig::default());
        let pf = opt.packing_fraction();
        assert!(pf > 0.0 && pf < 1.0);
    }
    #[test]
    fn test_void_fraction_empty_packing() {
        let analysis = VoidSpaceAnalysis::new(vec![], [1.0, 1.0, 1.0], 5);
        let vf = analysis.void_fraction_by_voxelization();
        assert!((vf - 1.0).abs() < 0.01);
    }
    #[test]
    fn test_void_fraction_full_packing() {
        let s = PackingSphere::new([0.5, 0.5, 0.5], 1.0, 0);
        let analysis = VoidSpaceAnalysis::new(vec![s], [1.0, 1.0, 1.0], 5);
        let vf = analysis.void_fraction_by_voxelization();
        assert!(vf < 0.5);
    }
    #[test]
    fn test_point_in_sphere_true() {
        let s = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let analysis = VoidSpaceAnalysis::new(vec![s], [2.0, 2.0, 2.0], 4);
        assert!(analysis.point_in_sphere([0.0, 0.0, 0.0]));
        assert!(analysis.point_in_sphere([0.9, 0.0, 0.0]));
    }
    #[test]
    fn test_point_in_sphere_false() {
        let s = PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0);
        let analysis = VoidSpaceAnalysis::new(vec![s], [3.0, 3.0, 3.0], 4);
        assert!(!analysis.point_in_sphere([2.0, 0.0, 0.0]));
    }
    #[test]
    fn test_pore_radius_in_void() {
        let s = PackingSphere::new([2.0, 2.0, 2.0], 0.5, 0);
        let analysis = VoidSpaceAnalysis::new(vec![s], [4.0, 4.0, 4.0], 4);
        let pr = analysis.pore_radius_at([3.0, 2.0, 2.0]);
        assert!(
            pr > 0.0,
            "pore radius at interior point should be positive, got {}",
            pr
        );
    }
    #[test]
    fn test_percolation_threshold() {
        let analysis = VoidSpaceAnalysis::new(vec![], [1.0, 1.0, 1.0], 4);
        assert!((analysis.estimate_percolation_threshold() - 0.3116).abs() < 0.001);
    }
    #[test]
    fn test_is_percolating_above() {
        let analysis = VoidSpaceAnalysis::new(vec![], [1.0, 1.0, 1.0], 4);
        assert!(analysis.is_percolating(0.5));
    }
    #[test]
    fn test_is_percolating_below() {
        let analysis = VoidSpaceAnalysis::new(vec![], [1.0, 1.0, 1.0], 4);
        assert!(!analysis.is_percolating(0.1));
    }
    #[test]
    fn test_specific_surface_area() {
        let s = PackingSphere::new([0.5, 0.5, 0.5], 0.1, 0);
        let analysis = VoidSpaceAnalysis::new(vec![s], [1.0, 1.0, 1.0], 4);
        let ssa = analysis.specific_surface_area();
        assert!(ssa > 0.0);
    }
    #[test]
    fn test_confine_box_inside() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Box,
            dimensions: [1.0, 1.0, 1.0],
            sphere_radius: 0.1,
            max_attempts: 100,
        };
        let conf = ConfiningGeometry::new(config);
        assert!(conf.point_inside_container([0.0, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_confine_box_outside() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Box,
            dimensions: [1.0, 1.0, 1.0],
            sphere_radius: 0.1,
            max_attempts: 100,
        };
        let conf = ConfiningGeometry::new(config);
        assert!(!conf.point_inside_container([0.95, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_confine_cylinder_inside() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Cylinder,
            dimensions: [1.0, 1.0, 0.0],
            sphere_radius: 0.1,
            max_attempts: 100,
        };
        let conf = ConfiningGeometry::new(config);
        assert!(conf.point_inside_container([0.0, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_confine_sphere_container_inside() {
        let config = ConfinementConfig {
            shape: ConfinementShape::SphericalContainer,
            dimensions: [2.0, 0.0, 0.0],
            sphere_radius: 0.1,
            max_attempts: 100,
        };
        let conf = ConfiningGeometry::new(config);
        assert!(conf.point_inside_container([0.0, 0.0, 0.0], 0.1));
        assert!(!conf.point_inside_container([1.95, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_confine_box_volume() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Box,
            dimensions: [1.0, 2.0, 3.0],
            sphere_radius: 0.1,
            max_attempts: 100,
        };
        let conf = ConfiningGeometry::new(config);
        assert!((conf.container_volume() - 48.0).abs() < 1e-10);
    }
    #[test]
    fn test_confine_cylinder_volume() {
        let r = 1.0;
        let h = 2.0;
        let config = ConfinementConfig {
            shape: ConfinementShape::Cylinder,
            dimensions: [r, h, 0.0],
            sphere_radius: 0.1,
            max_attempts: 100,
        };
        let conf = ConfiningGeometry::new(config);
        let expected = std::f64::consts::PI * r * r * 2.0 * h;
        assert!((conf.container_volume() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_confine_packing_run() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Box,
            dimensions: [2.0, 2.0, 2.0],
            sphere_radius: 0.15,
            max_attempts: 200,
        };
        let mut conf = ConfiningGeometry::new(config);
        let n = conf.run();
        assert!(n > 0);
        assert_eq!(conf.spheres.len(), n);
    }
    #[test]
    fn test_confine_no_overlap() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Box,
            dimensions: [2.0, 2.0, 2.0],
            sphere_radius: 0.2,
            max_attempts: 200,
        };
        let mut conf = ConfiningGeometry::new(config);
        conf.run();
        assert!(!has_overlap(&conf.spheres));
    }
    #[test]
    fn test_confine_packing_fraction_positive() {
        let config = ConfinementConfig {
            shape: ConfinementShape::Box,
            dimensions: [2.0, 2.0, 2.0],
            sphere_radius: 0.2,
            max_attempts: 200,
        };
        let mut conf = ConfiningGeometry::new(config);
        conf.run();
        assert!(conf.packing_fraction() > 0.0);
    }
    #[test]
    fn test_radial_distribution_function() {
        let mut packing = OrderedPacking::new(LatticeType::SimpleCubic, 1.0, [3, 3, 3]);
        packing.generate();
        let conf = ConfiningGeometry {
            config: ConfinementConfig::default(),
            spheres: packing.spheres.clone(),
            n_rejected: 0,
        };
        let (r, g) = conf.radial_distribution_function(20, 5.0);
        assert_eq!(r.len(), 20);
        assert_eq!(g.len(), 20);
    }
    #[test]
    fn test_compute_packing_fraction() {
        let spheres = vec![PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0)];
        let v = 8.0_f64.powi(3);
        let pf = compute_packing_fraction(&spheres, v);
        assert!(pf > 0.0 && pf < 1.0);
    }
    #[test]
    fn test_has_overlap_true() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([1.0, 0.0, 0.0], 1.0, 1),
        ];
        assert!(has_overlap(&spheres));
    }
    #[test]
    fn test_has_overlap_false() {
        let spheres = vec![
            PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0),
            PackingSphere::new([5.0, 0.0, 0.0], 1.0, 1),
        ];
        assert!(!has_overlap(&spheres));
    }
    #[test]
    fn test_structure_factor_single_sphere() {
        let spheres = vec![PackingSphere::new([0.0, 0.0, 0.0], 1.0, 0)];
        let q = vec![1.0, 2.0, 3.0];
        let sq = structure_factor(&spheres, &q);
        for s in &sq {
            assert!((s - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_structure_factor_empty() {
        let spheres: Vec<PackingSphere> = Vec::new();
        let q = vec![1.0];
        let sq = structure_factor(&spheres, &q);
        assert_eq!(sq, vec![0.0]);
    }
    #[test]
    fn test_coordination_number_distribution() {
        let mut packing = OrderedPacking::new(LatticeType::SimpleCubic, 1.0, [3, 3, 3]);
        packing.generate();
        let dist = coordination_number_distribution(&packing.spheres, 0.02);
        let total: usize = dist.iter().map(|(_, c)| c).sum();
        assert_eq!(total, packing.spheres.len());
    }
    #[test]
    fn test_polydisperse_gaussian_radius_bounds() {
        let params = SizeDistributionParams {
            mean_radius: 0.05,
            std_radius: 0.01,
            min_radius: 0.02,
            max_radius: 0.08,
            ..Default::default()
        };
        let packing = PolyDispersePacking::new(
            SizeDistribution::Gaussian,
            params.clone(),
            [2.0, 2.0, 2.0],
            1000,
        );
        let mut rng = rand::rng();
        for _ in 0..100 {
            let r = packing.sample_radius(&mut rng);
            assert!(r >= params.min_radius && r <= params.max_radius);
        }
    }
    #[test]
    fn test_polydisperse_powerlaw_radius_bounds() {
        let params = SizeDistributionParams {
            min_radius: 0.01,
            max_radius: 0.1,
            power_law_exponent: -3.0,
            ..Default::default()
        };
        let packing = PolyDispersePacking::new(
            SizeDistribution::PowerLaw,
            params.clone(),
            [2.0, 2.0, 2.0],
            1000,
        );
        let mut rng = rand::rng();
        for _ in 0..100 {
            let r = packing.sample_radius(&mut rng);
            assert!(
                r >= params.min_radius && r <= params.max_radius,
                "r={} out of bounds",
                r
            );
        }
    }
}
