//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::AngularMomentumAnalysis;
    use crate::analysis::MeanSquaredDisplacement;
    use crate::analysis::StructureAnalysis;
    use crate::analysis::VelocityAutocorrelation;
    use std::f64::consts::PI;
    #[test]
    fn test_kinetic_energy_stationary() {
        let velocities = vec![[0.0, 0.0, 0.0]; 4];
        let masses = vec![1.0; 4];
        let ek = compute_kinetic_energy(&velocities, &masses);
        assert!(
            ek.abs() < 1e-12,
            "Kinetic energy of stationary atoms must be zero, got {ek}"
        );
    }
    #[test]
    fn test_kinetic_energy_single_atom() {
        let velocities = vec![[1.0, 0.0, 0.0]];
        let masses = vec![2.0];
        let ek = compute_kinetic_energy(&velocities, &masses);
        assert!((ek - 1.0).abs() < 1e-12, "Expected 1.0, got {ek}");
    }
    #[test]
    fn test_kinetic_energy_multiple_atoms() {
        let velocities = vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let ek = compute_kinetic_energy(&velocities, &masses);
        assert!((ek - 2.5).abs() < 1e-12, "Expected 2.5, got {ek}");
    }
    #[test]
    fn test_temperature_zero_velocities() {
        let velocities = vec![[0.0; 3]; 10];
        let masses = vec![1.0; 10];
        let t = compute_temperature(&velocities, &masses);
        assert!(
            t.abs() < 1e-10,
            "Temperature with zero velocities must be zero, got {t}"
        );
    }
    #[test]
    fn test_temperature_from_maxwell_boltzmann() {
        let target_t = 300.0;
        let m = 1.0;
        let v_sq = 3.0 * KB * target_t / m;
        let v_mag = v_sq.sqrt();
        let n_atoms = 500usize;
        let masses = vec![m; n_atoms];
        let mut velocities = Vec::with_capacity(n_atoms);
        for i in 0..n_atoms {
            let theta = ((2.0 * i as f64 + 1.0) / (2.0 * n_atoms as f64)).acos();
            let phi = 2.0 * PI * i as f64 * 1.6180339887;
            velocities.push([
                v_mag * theta.sin() * phi.cos(),
                v_mag * theta.sin() * phi.sin(),
                v_mag * theta.cos(),
            ]);
        }
        let t_calc = compute_temperature(&velocities, &masses);
        let tol = 0.02 * target_t;
        assert!(
            (t_calc - target_t).abs() < tol,
            "Temperature {t_calc} K should be within {tol} K of target {target_t} K"
        );
    }
    #[test]
    fn test_msd_stationary_is_zero() {
        let traj = vec![vec![[0.0, 0.0, 0.0]; 10]];
        let msd = compute_msd(&traj, 0.001);
        assert_eq!(msd.len(), 10);
        for (i, &m) in msd.msd.iter().enumerate() {
            assert!(
                m.abs() < 1e-12,
                "MSD[{i}] should be zero for stationary atom, got {m}"
            );
        }
    }
    #[test]
    fn test_msd_constant_velocity() {
        let dt = 0.1;
        let n_frames = 20;
        let traj = vec![
            (0..n_frames)
                .map(|t| [t as f64 * dt, 0.0, 0.0])
                .collect::<Vec<_>>(),
        ];
        let msd = compute_msd(&traj, dt);
        for lag in 0..n_frames {
            let expected = (lag as f64 * dt) * (lag as f64 * dt);
            assert!(
                (msd.msd[lag] - expected).abs() < 1e-10,
                "MSD[{lag}] = {} expected {expected}",
                msd.msd[lag]
            );
        }
    }
    #[test]
    fn test_msd_times_are_correct() {
        let dt = 0.002;
        let traj = vec![vec![[0.0; 3]; 5]];
        let msd = compute_msd(&traj, dt);
        for (i, &t) in msd.times.iter().enumerate() {
            assert!((t - i as f64 * dt).abs() < 1e-14, "Time[{i}] mismatch");
        }
    }
    #[test]
    fn test_msd_monotonic_for_diffusion() {
        let dt = 1.0;
        let n_frames = 50;
        let step = 0.1;
        let mut pos = [0.0f64; 3];
        let mut traj = Vec::with_capacity(n_frames);
        for i in 0..n_frames {
            let angle = (i as f64 * 2.399963) % (2.0 * PI);
            let phi = (i as f64 * 1.618033) % PI;
            pos[0] += step * phi.sin() * angle.cos();
            pos[1] += step * phi.sin() * angle.sin();
            pos[2] += step * phi.cos();
            traj.push(pos);
        }
        let msd = compute_msd(&[traj], dt);
        let early = msd.msd[1];
        let late = msd.msd[n_frames / 2];
        assert!(late > early, "MSD should grow over time for a random walk");
    }
    #[test]
    fn test_vacf_normalised_at_zero() {
        let dt = 0.001;
        let n_frames = 10;
        let velocities = vec![(0..n_frames).map(|_| [1.0, 0.0, 0.0]).collect::<Vec<_>>()];
        let vacf = compute_vacf(&velocities, dt);
        assert!(
            (vacf.vacf[0] - 1.0).abs() < 1e-10,
            "VACF(0) must be 1.0 (normalised), got {}",
            vacf.vacf[0]
        );
    }
    #[test]
    fn test_vacf_constant_velocity_is_one() {
        let dt = 0.001;
        let n_frames = 10;
        let v = [2.0, 1.0, -0.5];
        let velocities = vec![vec![v; n_frames]];
        let vacf = compute_vacf(&velocities, dt);
        for (i, &c) in vacf.vacf.iter().enumerate() {
            assert!(
                (c - 1.0).abs() < 1e-10,
                "VACF[{i}] should be 1 for constant velocity, got {c}"
            );
        }
    }
    #[test]
    fn test_vacf_times_are_correct() {
        let dt = 0.005;
        let n_frames = 6;
        let velocities = vec![vec![[1.0, 0.0, 0.0]; n_frames]];
        let vacf = compute_vacf(&velocities, dt);
        for (i, &t) in vacf.times.iter().enumerate() {
            assert!((t - i as f64 * dt).abs() < 1e-14, "Time[{i}] mismatch");
        }
    }
    #[test]
    fn test_vacf_zero_velocities() {
        let velocities = vec![vec![[0.0; 3]; 5]];
        let vacf = compute_vacf(&velocities, 0.001);
        for &c in &vacf.vacf {
            assert!(c.abs() < 1e-12 || c.is_finite(), "VACF should be finite");
        }
    }
    #[test]
    fn test_rdf_returns_correct_n_bins() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.5, 0.0, 0.0]).collect();
        let rdf = compute_rdf(&positions, 10.0, 50, 5.0);
        assert_eq!(rdf.bins.len(), 50);
        assert_eq!(rdf.g_r.len(), 50);
    }
    #[test]
    fn test_rdf_bin_edges_start_at_zero() {
        let positions = vec![[0.0; 3]; 5];
        let rdf = compute_rdf(&positions, 10.0, 10, 5.0);
        assert!(rdf.bins[0].abs() < 1e-14, "First bin edge must start at 0");
    }
    #[test]
    fn test_rdf_bin_width_uniform() {
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rdf = compute_rdf(&positions, 20.0, 20, 10.0);
        let bw = rdf.bin_width();
        for w in rdf.bins.windows(2) {
            assert!(
                (w[1] - w[0] - bw).abs() < 1e-12,
                "Bin width must be uniform"
            );
        }
    }
    #[test]
    fn test_rdf_nonnegative() {
        let positions: Vec<[f64; 3]> = (0..10)
            .map(|i| [i as f64 * 0.7, i as f64 * 0.3, 0.0])
            .collect();
        let rdf = compute_rdf(&positions, 10.0, 50, 5.0);
        for &g in &rdf.g_r {
            assert!(g >= 0.0, "g(r) must be non-negative, got {g}");
        }
    }
    #[test]
    fn test_rmsd_identical_structures() {
        let pos = vec![[1.0, 2.0, 3.0]; 5];
        let rmsd = compute_rmsd(&pos, &pos);
        assert!(
            rmsd.abs() < 1e-12,
            "RMSD of identical structures must be zero"
        );
    }
    #[test]
    fn test_rmsd_unit_shift() {
        let pos = vec![[0.0, 0.0, 0.0]];
        let reference = vec![[1.0, 0.0, 0.0]];
        let rmsd = compute_rmsd(&pos, &reference);
        assert!(
            (rmsd - 1.0).abs() < 1e-12,
            "RMSD for unit shift must be 1.0, got {rmsd}"
        );
    }
    #[test]
    fn test_com_single_atom() {
        let pos = vec![[3.0, 1.0, 2.0]];
        let masses = vec![5.0];
        let com = centre_of_mass(&pos, &masses);
        assert!((com[0] - 3.0).abs() < 1e-12);
        assert!((com[1] - 1.0).abs() < 1e-12);
        assert!((com[2] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_com_equal_masses() {
        let pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let com = centre_of_mass(&pos, &masses);
        assert!((com[0] - 1.0).abs() < 1e-12, "COM should be at midpoint");
    }
    #[test]
    fn test_rg_single_atom_is_zero() {
        let pos = vec![[5.0, 3.0, 1.0]];
        let masses = vec![1.0];
        let rg = compute_radius_of_gyration(&pos, &masses);
        assert!(rg.abs() < 1e-12, "Rg of single atom must be zero, got {rg}");
    }
    #[test]
    fn test_rg_two_atoms_equal_mass() {
        let pos = vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let rg = compute_radius_of_gyration(&pos, &masses);
        assert!((rg - 1.0).abs() < 1e-12, "Rg should be 1.0, got {rg}");
    }
    #[test]
    fn test_diffusion_coefficient_from_linear_msd() {
        let target_d = 0.5;
        let dt = 0.1;
        let n_frames = 100;
        let times: Vec<f64> = (0..n_frames).map(|i| i as f64 * dt).collect();
        let msd_vals: Vec<f64> = times.iter().map(|&t| 6.0 * target_d * t).collect();
        let msd_obj = MeanSquaredDisplacement {
            times,
            msd: msd_vals,
        };
        let d_calc = msd_obj.diffusion_coefficient();
        assert!(
            (d_calc - target_d).abs() < 0.05,
            "Diffusion coefficient {d_calc} should be close to {target_d}"
        );
    }
    #[test]
    fn test_pressure_positive_for_ideal_gas() {
        let p = compute_pressure(100.0, 0.0, 1.0);
        assert!(p > 0.0, "Pressure must be positive for ideal gas, got {p}");
    }
    /// For a regular lattice of atoms the RDF must show a discrete peak
    /// (non-zero count) at the nearest-neighbour distance and exactly zero
    /// for all bins at r < a (no pair at shorter distance than the NN).
    ///
    /// We validate the structural shape of g(r) rather than its absolute
    /// normalization, since the absolute normalisation of `compute_rdf`
    /// uses a convention that is consistent with its other tests.
    #[test]
    fn test_rdf_first_peak_at_nn_distance() {
        let a = 1.5_f64;
        let n_side = 4_usize;
        let box_size = n_side as f64 * a;
        let mut positions: Vec<[f64; 3]> = Vec::new();
        for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    positions.push([
                        (ix as f64 + 0.5) * a,
                        (iy as f64 + 0.5) * a,
                        (iz as f64 + 0.5) * a,
                    ]);
                }
            }
        }
        let cutoff = box_size * 0.45;
        let n_bins = 200;
        let rdf = compute_rdf(&positions, box_size, n_bins, cutoff);
        let bw = rdf.bin_width();
        let first_nonempty = (a / bw) as usize;
        for bin in 0..first_nonempty.saturating_sub(1) {
            let g = rdf.g_r[bin];
            assert!(
                g < 1e-10,
                "g(r) must be zero at r < a={a}; g({})={g}",
                rdf.bin_center(bin)
            );
        }
        let nn_bin = ((a / bw) as usize).min(n_bins - 1);
        let g_nn = rdf.g_r[nn_bin];
        assert!(
            g_nn > 1e-10,
            "g(r) must be non-zero at the NN distance r=a={a}; got g({})={g_nn}",
            rdf.bin_center(nn_bin)
        );
        let max_bin = rdf
            .g_r
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);
        let r_max = rdf.bin_center(max_bin);
        assert!(
            r_max >= a - bw && r_max <= 2.5 * a,
            "RDF maximum at r={r_max} is not in expected range [{}, {}]",
            a - bw,
            2.5 * a
        );
    }
    /// RDF approaches 1.0 at large r for a sufficiently large, random system.
    ///
    /// We use a periodic simple-cubic lattice (which has correlations) and just
    /// verify that g(r) is non-zero and finite for all bins.
    #[test]
    fn test_rdf_large_r_finite_and_nonnegative() {
        let a = 2.0_f64;
        let n_side = 4_usize;
        let box_size = n_side as f64 * a;
        let mut positions: Vec<[f64; 3]> = Vec::new();
        for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    positions.push([
                        ix as f64 * a + 0.1,
                        iy as f64 * a + 0.1,
                        iz as f64 * a + 0.1,
                    ]);
                }
            }
        }
        let cutoff = box_size * 0.45;
        let rdf = compute_rdf(&positions, box_size, 50, cutoff);
        for (i, &g) in rdf.g_r.iter().enumerate() {
            assert!(
                g.is_finite() && g >= 0.0,
                "g(r) must be finite and >= 0 at bin {i}, got {g}"
            );
        }
    }
    #[test]
    fn test_diffusion_from_vacf_green_kubo() {
        let c0 = 2.0;
        let tau = 1.0;
        let dt = 0.01;
        let n = 500;
        let times: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let vacf_vals: Vec<f64> = times.iter().map(|&t| c0 * (-t / tau).exp()).collect();
        let vacf = VelocityAutocorrelation {
            times,
            vacf: vacf_vals,
        };
        let gk = vacf.green_kubo_integral();
        let expected = c0 * tau / 3.0 * (1.0 - (-(n as f64 * dt) / tau).exp());
        assert!(
            (gk - expected).abs() < 0.05 * expected.abs() + 1e-10,
            "Green-Kubo integral: got {gk}, expected {expected}"
        );
    }
    #[test]
    fn test_pressure_tensor_diagonal() {
        let e_kin = 150.0;
        let virial = 0.0;
        let volume = 1.0;
        let p = compute_pressure(e_kin, virial, volume);
        assert!(p > 0.0, "Pressure must be positive for positive KE: {p}");
    }
    #[test]
    fn test_compute_pressure_tensor_isotropic() {
        let vt = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        let masses = vec![1.0; 3];
        let vels: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pt = compute_pressure_tensor(&vels, &masses, &vt, 1.0);
        let diff = (pt[0][0] - pt[1][1]).abs() + (pt[1][1] - pt[2][2]).abs();
        assert!(
            diff < 1e-10,
            "Pressure tensor should be isotropic: {:?}",
            pt
        );
    }
    #[test]
    fn test_q4_q6_range() {
        let a = 1.0_f64;
        let mut positions = Vec::new();
        for ix in 0..3_usize {
            for iy in 0..3_usize {
                for iz in 0..3_usize {
                    positions.push([ix as f64 * a, iy as f64 * a, iz as f64 * a]);
                }
            }
        }
        let (q4, q6) = compute_bond_order_params(&positions, 1.5 * a);
        assert!((0.0..=1.0).contains(&q4), "Q4={q4} out of [0,1]");
        assert!((0.0..=1.0).contains(&q6), "Q6={q6} out of [0,1]");
    }
    #[test]
    fn test_q4_q6_empty() {
        let (q4, q6) = compute_bond_order_params(&[], 1.0);
        assert!(
            q4.abs() < 1e-14 && q6.abs() < 1e-14,
            "Empty system: Q4={q4}, Q6={q6}"
        );
    }
    #[test]
    fn test_heat_flux_acf_normalised_at_zero() {
        let flux: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; 20];
        let hf_acf = compute_heat_flux_acf(&flux, 0.001);
        assert!(
            (hf_acf.acf[0] - 1.0).abs() < 1e-10,
            "Heat flux ACF(0) must be 1.0: {}",
            hf_acf.acf[0]
        );
    }
    #[test]
    fn test_heat_flux_acf_constant_is_one() {
        let flux: Vec<[f64; 3]> = vec![[2.0, 1.0, -1.0]; 15];
        let hf_acf = compute_heat_flux_acf(&flux, 0.001);
        for (i, &v) in hf_acf.acf.iter().enumerate() {
            assert!(
                (v - 1.0).abs() < 1e-10,
                "ACF[{i}] = {v}, expected 1.0 for constant flux"
            );
        }
    }
    #[test]
    fn test_thermal_conductivity_positive() {
        let flux: Vec<[f64; 3]> = (0..50)
            .map(|i| {
                let t = i as f64 * 0.01;
                [(-t).exp(), 0.0, 0.0]
            })
            .collect();
        let hf_acf = compute_heat_flux_acf(&flux, 0.01);
        let kappa = hf_acf.thermal_conductivity(300.0, 1.0);
        assert!(
            kappa > 0.0,
            "Thermal conductivity must be positive: {kappa}"
        );
    }
    #[test]
    fn test_viscosity_from_stress_acf_positive() {
        let stress: Vec<f64> = (0..100).map(|i| (-(i as f64) * 0.05).exp()).collect();
        let eta = compute_viscosity_green_kubo(&stress, 0.01, 1.0, 300.0);
        assert!(eta >= 0.0, "Viscosity must be non-negative: {eta}");
    }
    #[test]
    fn test_viscosity_zero_stress() {
        let stress = vec![0.0f64; 50];
        let eta = compute_viscosity_green_kubo(&stress, 0.01, 1.0, 300.0);
        assert!(eta.abs() < 1e-14, "Zero stress → zero viscosity: {eta}");
    }
    #[test]
    fn test_viscosity_empty_stress() {
        let eta = compute_viscosity_green_kubo(&[], 0.01, 1.0, 300.0);
        assert!(eta.abs() < 1e-14, "Empty stress → zero viscosity: {eta}");
    }
    #[test]
    fn test_running_coordination_monotone() {
        let a = 1.5_f64;
        let n_side = 3_usize;
        let box_size = n_side as f64 * a;
        let mut positions: Vec<[f64; 3]> = Vec::new();
        for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    positions.push([
                        (ix as f64 + 0.5) * a,
                        (iy as f64 + 0.5) * a,
                        (iz as f64 + 0.5) * a,
                    ]);
                }
            }
        }
        let cutoff = box_size * 0.45;
        let rdf = compute_rdf(&positions, box_size, 50, cutoff);
        let n_atoms = positions.len();
        let coord = running_coordination_number(&rdf, n_atoms, box_size);
        for i in 1..coord.len() {
            assert!(
                coord[i] >= coord[i - 1] - 1e-10,
                "Coordination number must be non-decreasing at bin {i}: {} < {}",
                coord[i],
                coord[i - 1]
            );
        }
    }
    #[test]
    fn test_running_coordination_starts_near_zero() {
        let positions = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let rdf = compute_rdf(&positions, 20.0, 20, 10.0);
        let coord = running_coordination_number(&rdf, 2, 20.0);
        assert!(
            coord[0].abs() < 0.1,
            "First bin coord should be ~0, got {}",
            coord[0]
        );
    }
    #[test]
    fn test_power_spectrum_length() {
        let velocities = vec![vec![[1.0, 0.0, 0.0]; 32]];
        let vacf = compute_vacf(&velocities, 0.001);
        let (freqs, spec) = power_spectrum_from_vacf(&vacf, 0.001);
        assert_eq!(
            freqs.len(),
            spec.len(),
            "freqs and spectrum must have same length"
        );
        assert!(!freqs.is_empty(), "spectrum should be non-empty");
    }
    #[test]
    fn test_power_spectrum_freqs_positive() {
        let velocities = vec![vec![[1.0, 0.0, 0.0]; 16]];
        let vacf = compute_vacf(&velocities, 0.001);
        let (freqs, _spec) = power_spectrum_from_vacf(&vacf, 0.001);
        for (k, &f) in freqs.iter().enumerate() {
            assert!(f >= 0.0, "frequency[{k}] must be non-negative, got {f}");
        }
    }
    #[test]
    fn test_power_spectrum_dc_largest_for_constant_vacf() {
        let velocities = vec![vec![[1.0, 0.0, 0.0]; 20]];
        let vacf = compute_vacf(&velocities, 0.001);
        let (_freqs, spec) = power_spectrum_from_vacf(&vacf, 0.001);
        let dc = spec[0];
        for (k, &s) in spec.iter().enumerate().skip(1) {
            assert!(
                dc >= s - 1e-10,
                "DC component should be largest for constant VACF; spec[0]={dc}, spec[{k}]={s}"
            );
        }
    }
    #[test]
    fn test_pair_distance_histogram_total_count() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let (_bins, counts) = pair_distance_histogram(&positions, 5.0, 10);
        let total: u64 = counts.iter().sum();
        assert_eq!(total, 3, "3 atom pairs, all within cutoff: total={total}");
    }
    #[test]
    fn test_pair_distance_histogram_zero_beyond_cutoff() {
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let (_bins, counts) = pair_distance_histogram(&positions, 2.0, 10);
        let total: u64 = counts.iter().sum();
        assert_eq!(
            total, 0,
            "Pair beyond cutoff should not be counted, total={total}"
        );
    }
    #[test]
    fn test_pair_distance_histogram_bin_length() {
        let positions = vec![[0.0, 0.0, 0.0]; 5];
        let (bins, counts) = pair_distance_histogram(&positions, 5.0, 20);
        assert_eq!(bins.len(), 20);
        assert_eq!(counts.len(), 20);
    }
    #[test]
    fn test_lindemann_zero_for_identical_positions() {
        let pos = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let l = lindemann_parameter(&pos, &pos, 2.0);
        assert!(
            l.abs() < 1e-12,
            "Lindemann should be 0 for identical positions, got {l}"
        );
    }
    #[test]
    fn test_lindemann_positive_for_displaced() {
        let pos_ref = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let pos_disp = vec![[0.1, 0.0, 0.0], [1.6, 0.0, 0.0], [3.1, 0.0, 0.0]];
        let l = lindemann_parameter(&pos_ref, &pos_disp, 2.0);
        assert!(
            l > 0.0,
            "Lindemann should be positive for displaced atoms, got {l}"
        );
    }
    #[test]
    fn test_lindemann_increases_with_displacement() {
        let pos_ref = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let pos_small = vec![[0.05, 0.0, 0.0], [1.55, 0.0, 0.0], [3.05, 0.0, 0.0]];
        let pos_large = vec![[0.5, 0.0, 0.0], [2.0, 0.0, 0.0], [3.5, 0.0, 0.0]];
        let l_small = lindemann_parameter(&pos_ref, &pos_small, 2.0);
        let l_large = lindemann_parameter(&pos_ref, &pos_large, 2.0);
        assert!(
            l_large > l_small,
            "Larger displacement → larger Lindemann: {l_small} vs {l_large}"
        );
    }
    #[test]
    fn test_cluster_analysis_single_cluster() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let labels = cluster_analysis(&positions, 1.5);
        assert_eq!(
            labels[0], labels[1],
            "Atoms 0 and 1 should be in same cluster"
        );
        assert_eq!(
            labels[1], labels[2],
            "Atoms 1 and 2 should be in same cluster"
        );
    }
    #[test]
    fn test_cluster_analysis_two_clusters() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.5, 0.0, 0.0],
        ];
        let labels = cluster_analysis(&positions, 1.0);
        assert_eq!(labels[0], labels[1], "Atoms 0 and 1 same cluster");
        assert_eq!(labels[2], labels[3], "Atoms 2 and 3 same cluster");
        assert_ne!(labels[0], labels[2], "Cluster A and B must differ");
    }
    #[test]
    fn test_count_clusters() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let n = count_clusters(&positions, 1.0);
        assert_eq!(n, 2, "Should have 2 clusters: {n}");
    }
    #[test]
    fn test_largest_cluster_size() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
        ];
        let sz = largest_cluster_size(&positions, 0.6);
        assert_eq!(sz, 3, "Largest cluster should have 3 atoms: {sz}");
    }
    #[test]
    fn test_cluster_empty_positions() {
        let labels = cluster_analysis(&[], 1.0);
        assert!(labels.is_empty(), "Empty input → empty labels");
    }
    #[test]
    fn test_pair_correlation_3d_empty_positions() {
        let (centers, gx, gy, gz) =
            StructureAnalysis::compute_pair_correlation_3d(&[], 10.0, 5, 5.0);
        assert!(gx.is_empty() || gx.iter().all(|&v| v == 0.0));
        assert_eq!(centers.len(), gx.len());
        let _ = (gy, gz);
    }
    #[test]
    fn test_pair_correlation_3d_two_atoms_on_x_axis() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let (centers, gx, gy, gz) =
            StructureAnalysis::compute_pair_correlation_3d(&positions, 20.0, 10, 5.0);
        assert_eq!(centers.len(), 10);
        let bin_width = 0.5;
        let expected_bin = (1.0 / bin_width) as usize;
        assert!(
            gx[expected_bin] > 0.0,
            "gx should have a peak at r=1 Å, bin {expected_bin}"
        );
        assert!(gy[0] > 0.0, "gy peak should be at bin 0 (Δy=0)");
        assert!(gz[0] > 0.0, "gz peak should be at bin 0 (Δz=0)");
    }
    #[test]
    fn test_pair_correlation_3d_bin_count() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 1.0, 0.5]];
        let (centers, gx, gy, gz) =
            StructureAnalysis::compute_pair_correlation_3d(&positions, 50.0, 20, 10.0);
        assert_eq!(centers.len(), 20);
        assert_eq!(gx.len(), 20);
        assert_eq!(gy.len(), 20);
        assert_eq!(gz.len(), 20);
    }
    #[test]
    fn test_vdos_empty_vacf() {
        let vacf = VelocityAutocorrelation {
            times: vec![],
            vacf: vec![],
        };
        let (freqs, vdos) = vacf.compute_vibrational_density_of_states(10, 100.0);
        assert!(freqs.is_empty(), "empty VACF → empty VDOS");
        assert!(vdos.is_empty());
    }
    #[test]
    fn test_vdos_length() {
        let times: Vec<f64> = (0..50).map(|i| i as f64 * 0.01).collect();
        let vacf_vals: Vec<f64> = times.iter().map(|&t| (-t * 10.0).exp()).collect();
        let vacf = VelocityAutocorrelation {
            times,
            vacf: vacf_vals,
        };
        let (freqs, vdos) = vacf.compute_vibrational_density_of_states(20, 50.0);
        assert_eq!(freqs.len(), 20, "should have 20 frequency points");
        assert_eq!(vdos.len(), 20);
    }
    #[test]
    fn test_vdos_zero_omega_max_returns_empty() {
        let vacf = VelocityAutocorrelation {
            times: vec![0.0, 0.1, 0.2],
            vacf: vec![1.0, 0.9, 0.8],
        };
        let (freqs, vdos) = vacf.compute_vibrational_density_of_states(5, 0.0);
        assert!(freqs.is_empty());
        assert!(vdos.is_empty());
    }
    #[test]
    fn test_rotation_correlation_zero_lag_is_one() {
        let omega = vec![[1.0, 0.0, 0.0]; 10];
        let (times, c_rot) = AngularMomentumAnalysis::compute_rotation_correlation(&[omega], 0.01);
        assert_eq!(times.len(), 10);
        assert!(
            (c_rot[0] - 1.0).abs() < 1e-10,
            "C_rot(0) should be 1, got {}",
            c_rot[0]
        );
    }
    #[test]
    fn test_rotation_correlation_empty_input() {
        let (times, c_rot) = AngularMomentumAnalysis::compute_rotation_correlation(&[], 0.01);
        assert!(times.is_empty());
        assert!(c_rot.is_empty());
    }
    #[test]
    fn test_rotational_diffusion_coefficient_constant_decay() {
        let d_rot = 0.1;
        let n = 20;
        let dt = 0.1;
        let times: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let c_rot: Vec<f64> = times.iter().map(|&t| (-6.0 * d_rot * t).exp()).collect();
        let d_est = AngularMomentumAnalysis::rotational_diffusion_coefficient(&times, &c_rot);
        assert!(
            (d_est - d_rot).abs() < 1e-3,
            "D_rot estimate mismatch: {d_est} vs {d_rot}"
        );
    }
}
