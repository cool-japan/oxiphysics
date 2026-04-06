//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Boltzmann constant for sampling module (J/K).
pub(super) const KB_SAMP: f64 = 1.380_649e-23;
#[cfg(test)]
mod tests {
    use super::super::types::*;

    use crate::enhanced_sampling::ReplicaExchange;
    use crate::sampling::AnnealingSchedule;
    use crate::sampling::SimulatedAnnealing;
    use crate::sampling::TransitionPath;
    #[test]
    fn test_umbrella_bias_force_zero_at_center() {
        let umb = UmbrellaSampling::new(2.0, 100.0);
        let force = umb.bias_force(2.0);
        assert!(
            force.abs() < 1e-12,
            "Force at center must be zero, got {force}"
        );
    }
    #[test]
    fn test_umbrella_bias_force_restoring() {
        let umb = UmbrellaSampling::new(1.0, 50.0);
        let force = umb.bias_force(1.5);
        assert!(
            force < 0.0,
            "Force must be negative (restoring) when xi > xi0, got {force}"
        );
    }
    #[test]
    fn test_umbrella_collect_statistics() {
        let mut umb = UmbrellaSampling::new(3.0, 1.0);
        for x in [1.0, 2.0, 3.0, 4.0, 5.0] {
            umb.collect(x);
        }
        let (mean, variance) = umb.statistics();
        assert!((mean - 3.0).abs() < 1e-10, "Mean should be 3, got {mean}");
        assert!(
            (variance - 2.0).abs() < 1e-10,
            "Variance should be 2, got {variance}"
        );
    }
    #[test]
    fn test_umbrella_histogram_bins() {
        let mut umb = UmbrellaSampling::new(0.5, 10.0);
        for i in 0..100 {
            umb.collect(i as f64 / 99.0);
        }
        let n_bins = 10;
        let (centers, energies) = umb.free_energy_histogram(n_bins, 1.0);
        assert_eq!(centers.len(), n_bins, "Should have {n_bins} bin centers");
        assert_eq!(
            energies.len(),
            n_bins,
            "Should have {n_bins} free energy values"
        );
    }
    #[test]
    fn test_umbrella_free_energy_shape() {
        let xi0 = 0.0;
        let mut umb = UmbrellaSampling::new(xi0, 1.0);
        let n = 2000usize;
        for i in 0..n {
            let sum: f64 = (0..12)
                .map(|j| {
                    let t = (i * 12 + j) as f64;
                    ((t * 1664525.0 + 1013904223.0).sin().abs() % 1.0).abs()
                })
                .sum::<f64>()
                - 6.0;
            umb.collect(sum * 0.2);
        }
        let (centers, energies) = umb.free_energy_histogram(20, 1.0);
        assert!(!centers.is_empty(), "Should have histogram bins");
        let min_idx = energies
            .iter()
            .enumerate()
            .filter(|&(_, &e)| e.is_finite())
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .expect("Should find minimum free energy bin");
        let min_xi = centers[min_idx];
        assert!(
            (min_xi - xi0).abs() < 0.5,
            "Free energy minimum at {min_xi} should be near xi0={xi0}"
        );
    }
    #[test]
    fn test_wham_single_window() {
        let xi0 = 1.0;
        let k = 10.0;
        let k_t = 1.0;
        let mut umb = UmbrellaSampling::new(xi0, k);
        for i in 0..200 {
            umb.collect(xi0 + (i as f64 - 100.0) * 0.01);
        }
        let (hcenters, _henergies) = umb.free_energy_histogram(10, k_t);
        let mut umb2 = UmbrellaSampling::new(xi0, k);
        for i in 0..200 {
            umb2.collect(xi0 + (i as f64 - 100.0) * 0.01);
        }
        let mut wham = WhamAnalysis::new(k_t);
        wham.add_window(umb2);
        let (wcenters, _wenergies) = wham.compute_pmf(10, 100);
        assert_eq!(hcenters.len(), 10, "Histogram should have 10 bins");
        assert_eq!(wcenters.len(), 10, "WHAM PMF should have 10 bins");
        for (hc, wc) in hcenters.iter().zip(wcenters.iter()) {
            assert!(
                (hc - wc).abs() < 1e-10,
                "Bin centers should match: histogram={hc}, wham={wc}"
            );
        }
    }
    #[test]
    fn test_replica_exchange_creation() {
        let re = ReplicaExchange::new(vec![300.0, 350.0, 400.0]);
        assert_eq!(re.n_replicas(), 3);
        assert_eq!(re.n_accepted, 0);
        assert_eq!(re.n_attempts, 0);
    }
    #[test]
    fn test_replica_exchange_favorable_swap() {
        let mut re = ReplicaExchange::new(vec![300.0, 400.0]);
        re.set_energy(0, 1e-19);
        re.set_energy(1, -1e-19);
        let accepted = re.try_swap(0, 1, 0.999);
        assert!(accepted, "Favorable swap should be accepted");
        assert_eq!(re.n_accepted, 1);
    }
    #[test]
    fn test_replica_exchange_acceptance_rate() {
        let mut re = ReplicaExchange::new(vec![300.0, 400.0]);
        re.set_energy(0, 1e-19);
        re.set_energy(1, -1e-19);
        for _ in 0..10 {
            re.try_swap(0, 1, 0.5);
        }
        let rate = re.acceptance_rate();
        assert!(
            (0.0..=1.0).contains(&rate),
            "Rate must be in [0,1], got {rate}"
        );
    }
    #[test]
    fn test_sa_linear_initial() {
        let sa = SimulatedAnnealing::new(1000.0, 100.0, 100, AnnealingSchedule::Linear);
        let t = sa.temperature();
        assert!(
            (t - 1000.0).abs() < 1e-10,
            "Initial temp should be 1000, got {t}"
        );
    }
    #[test]
    fn test_sa_linear_final() {
        let mut sa = SimulatedAnnealing::new(1000.0, 100.0, 100, AnnealingSchedule::Linear);
        for _ in 0..100 {
            sa.step();
        }
        let t = sa.temperature();
        assert!(
            (t - 100.0).abs() < 1e-10,
            "Final temp should be 100, got {t}"
        );
    }
    #[test]
    fn test_sa_exponential_monotone() {
        let mut sa = SimulatedAnnealing::new(1000.0, 10.0, 50, AnnealingSchedule::Exponential);
        let mut prev = sa.temperature();
        for _ in 0..50 {
            sa.step();
            let t = sa.temperature();
            assert!(
                t <= prev + 1e-10,
                "Temperature should decrease monotonically"
            );
            prev = t;
        }
    }
    #[test]
    fn test_sa_metropolis_accept_downhill() {
        let sa = SimulatedAnnealing::new(1000.0, 100.0, 100, AnnealingSchedule::Linear);
        assert!(sa.metropolis_accept(-1.0, 0.99));
    }
    #[test]
    fn test_sa_reset() {
        let mut sa = SimulatedAnnealing::new(1000.0, 100.0, 100, AnnealingSchedule::Linear);
        for _ in 0..50 {
            sa.step();
        }
        sa.reset();
        assert_eq!(sa.current_step, 0);
        assert!((sa.temperature() - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_wl_creation() {
        let wl = WangLandau::new(-10.0, 10.0, 20, 0.8);
        assert_eq!(wl.n_bins(), 20);
        assert_eq!(wl.bin_edges.len(), 21);
        assert!((wl.ln_f - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_wl_bin_index() {
        let wl = WangLandau::new(0.0, 10.0, 10, 0.8);
        assert_eq!(wl.bin_index(0.5), Some(0));
        assert_eq!(wl.bin_index(9.5), Some(9));
        assert_eq!(wl.bin_index(-0.1), None);
        assert_eq!(wl.bin_index(10.0), None);
    }
    #[test]
    fn test_wl_visit_updates_ln_g() {
        let mut wl = WangLandau::new(0.0, 10.0, 10, 0.8);
        wl.visit(0.5);
        assert!((wl.ln_g[0] - 1.0).abs() < 1e-15);
        assert_eq!(wl.histogram[0], 1);
        wl.visit(0.5);
        assert!((wl.ln_g[0] - 2.0).abs() < 1e-15);
        assert_eq!(wl.histogram[0], 2);
    }
    #[test]
    fn test_wl_flatness_check() {
        let mut wl = WangLandau::new(0.0, 10.0, 5, 0.8);
        for i in 0..5 {
            let e = i as f64 + 0.5;
            for _ in 0..10 {
                wl.visit(e * 2.0);
            }
        }
        assert!(wl.is_flat(), "Uniform histogram should be flat");
    }
    #[test]
    fn test_wl_update_modification_factor() {
        let mut wl = WangLandau::new(0.0, 10.0, 5, 0.8);
        wl.visit(0.5);
        let cont = wl.update_modification_factor();
        assert!(cont, "Should continue (ln_f still large)");
        assert!((wl.ln_f - 0.5).abs() < 1e-15);
        assert_eq!(wl.histogram[0], 0, "Histogram should be reset");
    }
    #[test]
    fn test_wl_bin_centers() {
        let wl = WangLandau::new(0.0, 10.0, 5, 0.8);
        let centers = wl.bin_centers();
        assert_eq!(centers.len(), 5);
        assert!((centers[0] - 1.0).abs() < 1e-12);
        assert!((centers[4] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_transition_path_creation() {
        let path = TransitionPath::new(vec![0.1, 0.5, 0.9]);
        assert_eq!(path.len(), 3);
        assert!(!path.is_empty());
    }
    #[test]
    fn test_transition_path_connects_basins() {
        let path = TransitionPath::new(vec![0.1, 0.4, 0.6, 0.9]);
        assert!(path.connects_basins(0.2, 0.8));
        let path2 = TransitionPath::new(vec![0.1, 0.3, 0.5]);
        assert!(!path2.connects_basins(0.2, 0.8));
    }
    #[test]
    fn test_transition_path_reverse_direction() {
        let path = TransitionPath::new(vec![0.9, 0.6, 0.1]);
        assert!(path.connects_basins(0.2, 0.8));
    }
    #[test]
    fn test_tps_add_path() {
        let mut tps = TransitionPathSampling::new(0.2, 0.8);
        let good_path = TransitionPath::new(vec![0.1, 0.5, 0.9]);
        assert!(tps.add_path(good_path));
        assert_eq!(tps.paths.len(), 1);
        let bad_path = TransitionPath::new(vec![0.1, 0.3, 0.5]);
        assert!(!tps.add_path(bad_path));
        assert_eq!(tps.paths.len(), 1);
    }
    #[test]
    fn test_tps_shooting_move_accept() {
        let mut tps = TransitionPathSampling::new(0.2, 0.8);
        let initial = TransitionPath::new(vec![0.1, 0.5, 0.9]);
        tps.add_path(initial);
        let new_path = TransitionPath::new(vec![0.15, 0.45, 0.85]);
        let accepted = tps.shooting_move(0, new_path);
        assert!(accepted);
        assert_eq!(tps.accepted_shoots, 1);
        assert_eq!(tps.attempted_shoots, 1);
    }
    #[test]
    fn test_tps_shooting_move_reject() {
        let mut tps = TransitionPathSampling::new(0.2, 0.8);
        let initial = TransitionPath::new(vec![0.1, 0.5, 0.9]);
        tps.add_path(initial);
        let bad_path = TransitionPath::new(vec![0.1, 0.3, 0.5]);
        let accepted = tps.shooting_move(0, bad_path);
        assert!(!accepted);
        assert_eq!(tps.accepted_shoots, 0);
        assert_eq!(tps.attempted_shoots, 1);
    }
    #[test]
    fn test_tps_acceptance_rate() {
        let mut tps = TransitionPathSampling::new(0.2, 0.8);
        let initial = TransitionPath::new(vec![0.1, 0.5, 0.9]);
        tps.add_path(initial);
        let good = TransitionPath::new(vec![0.15, 0.45, 0.85]);
        tps.shooting_move(0, good);
        let bad = TransitionPath::new(vec![0.1, 0.3]);
        tps.shooting_move(0, bad);
        let rate = tps.acceptance_rate();
        assert!((rate - 0.5).abs() < 1e-10, "Expected rate 0.5, got {rate}");
    }
}
#[cfg(test)]
mod tests_new_sampling {
    use super::super::types::*;

    #[test]
    fn test_string_method_creation() {
        let sm = StringMethod::new(vec![0.0, 0.0], vec![1.0, 1.0], 5, 0.01);
        assert_eq!(sm.images.len(), 5);
        assert!((sm.images[0].coords[0] - 0.0).abs() < 1e-12);
        assert!((sm.images[4].coords[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_string_method_linear_interpolation() {
        let sm = StringMethod::new(vec![0.0], vec![4.0], 5, 0.01);
        for (i, im) in sm.images.iter().enumerate() {
            let expected = i as f64;
            assert!(
                (im.coords[0] - expected).abs() < 1e-10,
                "image {i} at {} expected {}",
                im.coords[0],
                expected
            );
        }
    }
    #[test]
    fn test_string_reparameterise() {
        let mut sm = StringMethod::new(vec![0.0], vec![10.0], 6, 0.01);
        sm.images[2].coords[0] = 1.5;
        sm.reparameterise();
        for i in 1..sm.images.len() {
            assert!(
                sm.images[i].coords[0] >= sm.images[i - 1].coords[0] - 1e-10,
                "images should be monotonically ordered after reparameterisation"
            );
        }
    }
    #[test]
    fn test_string_path_length() {
        let sm = StringMethod::new(vec![0.0], vec![1.0], 3, 0.01);
        let length = sm.path_length();
        assert!(
            (length - 1.0).abs() < 1e-10,
            "path length should be 1.0, got {length}"
        );
    }
    #[test]
    fn test_string_step() {
        let mut sm = StringMethod::new(vec![0.0], vec![2.0], 4, 0.001);
        let force_fn = |coords: &[f64]| -> (Vec<f64>, f64) {
            let x = coords[0];
            (vec![-2.0 * x], x * x)
        };
        sm.step(force_fn);
        assert_eq!(sm.iterations, 1);
    }
    #[test]
    fn test_string_transition_state() {
        let mut sm = StringMethod::new(vec![0.0], vec![4.0], 5, 0.01);
        for (i, en) in [0.0, 1.0, 3.0, 2.0, 1.0].iter().enumerate() {
            sm.images[i].energy = *en;
        }
        let ts = sm.transition_state_index();
        assert_eq!(ts, 2, "transition state should be at index 2");
    }
    #[test]
    fn test_neb_creation() {
        let neb = NudgedElasticBand::new(vec![0.0, 0.0], vec![1.0, 1.0], 5, 1.0, 0.01, 1e-4, false);
        assert_eq!(neb.n_images(), 5);
        assert!((neb.images[0].coords[0] - 0.0).abs() < 1e-12);
        assert!((neb.images[4].coords[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_neb_step_returns_force() {
        let mut neb = NudgedElasticBand::new(vec![0.0], vec![2.0], 5, 1.0, 0.001, 1e-4, false);
        let max_f = neb.step(|coords| {
            let x = coords[0];
            (vec![-2.0 * x], x * x)
        });
        assert!(max_f >= 0.0, "max force should be non-negative");
        assert_eq!(neb.iterations, 1);
    }
    #[test]
    fn test_neb_activation_energy() {
        let mut neb = NudgedElasticBand::new(vec![0.0], vec![4.0], 5, 1.0, 0.001, 1e-4, false);
        neb.images[0].energy = 0.0;
        neb.images[1].energy = 1.0;
        neb.images[2].energy = 2.0;
        neb.images[3].energy = 1.5;
        neb.images[4].energy = 1.0;
        let ea = neb.activation_energy();
        assert!(
            (ea - 2.0).abs() < 1e-12,
            "activation energy should be 2.0, got {ea}"
        );
    }
    #[test]
    fn test_neb_convergence_check() {
        let neb = NudgedElasticBand::new(vec![0.0], vec![1.0], 3, 1.0, 0.001, 100.0, false);
        assert!(
            neb.is_converged(),
            "should be converged with zero forces and large tol"
        );
    }
    #[test]
    fn test_neb_climbing_image() {
        let mut neb = NudgedElasticBand::new(vec![0.0], vec![4.0], 5, 1.0, 0.001, 1e-4, true);
        let _max_f = neb.step(|coords| {
            let x = coords[0];
            let e = (x - 2.0).powi(2) * (-(x - 2.0).powi(2) / 4.0 + 1.0);
            let de_dx = 2.0 * (x - 2.0) * (1.0 - (x - 2.0).powi(2) / 2.0);
            (vec![-de_dx], e)
        });
        assert_eq!(neb.iterations, 1);
    }
    #[test]
    fn test_path_cv_endpoints() {
        let images = vec![vec![0.0, 0.0], vec![0.5, 0.5], vec![1.0, 1.0]];
        let pcv = PathCV::new(images, 10.0);
        let (s, _z) = pcv.compute(&[0.0, 0.0]);
        assert!(s < 0.1, "s at start should be ~0, got {s}");
        let (s2, _z2) = pcv.compute(&[1.0, 1.0]);
        assert!(s2 > 0.9, "s at end should be ~1, got {s2}");
    }
    #[test]
    fn test_path_cv_midpoint() {
        let images = vec![vec![0.0], vec![0.5], vec![1.0]];
        let pcv = PathCV::new(images, 10.0);
        let (s, z) = pcv.compute(&[0.5]);
        assert!(
            (s - 0.5).abs() < 0.1,
            "s at midpoint should be ~0.5, got {s}"
        );
        assert!(z.is_finite(), "z should be finite, got {z}");
    }
    #[test]
    fn test_path_cv_n_images() {
        let images: Vec<Vec<f64>> = (0..7).map(|i| vec![i as f64]).collect();
        let pcv = PathCV::new(images, 1.0);
        assert_eq!(pcv.n_images(), 7);
    }
    #[test]
    fn test_committor_add_point() {
        let mut ca = CommittorAnalysis::new(0.2, 0.8);
        ca.add_point(0.3, 0.1);
        ca.add_point(0.5, 0.5);
        ca.add_point(0.7, 0.9);
        assert_eq!(ca.n_points(), 3);
    }
    #[test]
    fn test_committor_transition_state() {
        let mut ca = CommittorAnalysis::new(0.2, 0.8);
        ca.add_point(0.3, 0.1);
        ca.add_point(0.5, 0.49);
        ca.add_point(0.6, 0.7);
        let ts = ca.transition_state_xi().expect("should find TS");
        assert!((ts - 0.5).abs() < 0.15, "TS should be near 0.5, got {ts}");
    }
    #[test]
    fn test_committor_from_trajectories() {
        let paths: Vec<(f64, bool)> = vec![
            (0.3, false),
            (0.3, false),
            (0.3, false),
            (0.3, true),
            (0.5, true),
            (0.5, false),
            (0.5, true),
            (0.5, false),
            (0.7, true),
            (0.7, true),
            (0.7, true),
            (0.7, false),
        ];
        let mut ca = CommittorAnalysis::new(0.2, 0.8);
        ca.from_trajectories(&paths, 5);
        assert!(
            ca.n_points() > 0,
            "should have data points from trajectories"
        );
    }
    #[test]
    fn test_committor_mean() {
        let mut ca = CommittorAnalysis::new(0.2, 0.8);
        ca.add_point(0.3, 0.0);
        ca.add_point(0.5, 0.5);
        ca.add_point(0.7, 1.0);
        let mean = ca.mean_committor();
        assert!(
            (mean - 0.5).abs() < 1e-10,
            "mean committor should be 0.5, got {mean}"
        );
    }
    #[test]
    fn test_committor_clamp() {
        let mut ca = CommittorAnalysis::new(0.0, 1.0);
        ca.add_point(0.5, 1.5);
        ca.add_point(0.5, -0.1);
        assert!((ca.data[0].1 - 1.0).abs() < 1e-10);
        assert!((ca.data[1].1 - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_committor_empty() {
        let ca = CommittorAnalysis::new(0.2, 0.8);
        assert_eq!(ca.n_points(), 0);
        assert!(ca.transition_state_xi().is_none());
        assert!((ca.mean_committor() - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_nested_sampling_empty_evidence() {
        let ns = NestedSampling::new(100);
        let lz = ns.compute_evidence();
        assert!(
            lz.is_infinite() && lz < 0.0,
            "empty NS should return -inf evidence"
        );
    }
    #[test]
    fn test_nested_sampling_monotone_log_volumes() {
        let mut ns = NestedSampling::new(50);
        for i in 0..5 {
            ns.record_dead_point(-(i as f64));
        }
        for w in ns.log_prior_volumes.windows(2) {
            assert!(
                w[0] > w[1],
                "log volumes should decrease: {} > {}",
                w[0],
                w[1]
            );
        }
    }
    #[test]
    fn test_nested_sampling_evidence_finite() {
        let mut ns = NestedSampling::new(20);
        for _ in 0..40 {
            ns.record_dead_point(0.0);
        }
        let lz = ns.compute_evidence();
        assert!(
            lz.is_finite(),
            "evidence should be finite for flat likelihood, got {lz}"
        );
        assert!(
            lz < 0.1,
            "log Z for unit likelihood should be close to 0, got {lz}"
        );
    }
    #[test]
    fn test_transition_matrix_symmetric_entropy_diff() {
        let mut tm = TransitionMatrix::new(3);
        for _ in 0..10 {
            tm.record_transition(0, 1);
            tm.record_transition(1, 0);
        }
        let ds = tm.compute_entropy_difference(0, 1).unwrap();
        assert!(ds.abs() < 1e-10, "symmetric counts → ΔS = 0, got {ds}");
    }
    #[test]
    fn test_transition_matrix_no_counts_returns_none() {
        let tm = TransitionMatrix::new(4);
        assert!(
            tm.compute_entropy_difference(0, 1).is_none(),
            "no counts should return None"
        );
    }
    #[test]
    fn test_transition_matrix_update_entropy() {
        let mut tm = TransitionMatrix::new(3);
        for _ in 0..20 {
            tm.record_transition(0, 1);
        }
        for _ in 0..10 {
            tm.record_transition(1, 0);
        }
        tm.update_entropy_estimates();
        assert!(
            (tm.ln_g[0] - 0.0).abs() < 1e-10,
            "ln_g[0] should be 0 reference"
        );
        assert!(
            (tm.ln_g[1] - 2_f64.ln()).abs() < 1e-10,
            "ln_g[1] should be ln 2, got {}",
            tm.ln_g[1]
        );
    }
    #[test]
    fn test_muca_bin_index_in_range() {
        let muca = MulticanonicalSimulation::new(0.0, 10.0, 10);
        assert_eq!(muca.bin_index(0.5), Some(0));
        assert_eq!(muca.bin_index(9.9), Some(9));
        assert!(muca.bin_index(-0.1).is_none());
        assert!(muca.bin_index(10.0).is_none());
    }
    #[test]
    fn test_muca_compute_dos_flat_histogram() {
        let mut muca = MulticanonicalSimulation::new(0.0, 5.0, 5);
        for i in 0..5 {
            for _ in 0..10 {
                let e = i as f64 + 0.5;
                muca.visit(e);
            }
        }
        muca.compute_density_of_states();
        for (i, &lng) in muca.ln_dos.iter().enumerate() {
            assert!(
                lng.abs() < 1e-10,
                "flat histogram → ln_dos[{i}] = 0, got {lng}"
            );
        }
    }
    #[test]
    fn test_muca_weight_uniform_dos() {
        let muca = MulticanonicalSimulation::new(0.0, 4.0, 4);
        assert!(
            (muca.weight(0.5) - 1.0).abs() < 1e-10,
            "uniform dos → weight = 1"
        );
        assert!(muca.weight(-1.0).abs() < 1e-10, "out-of-range → weight = 0");
    }
}
