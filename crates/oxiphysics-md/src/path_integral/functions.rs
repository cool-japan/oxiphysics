//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

/// Reduced Planck constant (J*s).
pub(super) const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant (J K^-1).
pub(super) const KB: f64 = 1.380_649e-23;
/// Build an orthogonal P x P normal mode transform matrix.
///
/// Row 0: centroid mode (1/sqrt(P) for all n)
/// For P even:
///   Row P/2: alternating mode ((-1)^n / sqrt(P))
///   Rows k (1 <= k < P/2): sqrt(2/P) * cos(2*pi*k*n/P)
///   Rows P-k: sqrt(2/P) * sin(2*pi*k*n/P)
/// For P odd:
///   Rows k (1 <= k <= (P-1)/2): sqrt(2/P) * cos(2*pi*k*n/P)
///   Rows P-k: sqrt(2/P) * sin(2*pi*k*n/P)
///
/// Returns a flat vector of length P*P stored row-major.
pub(super) fn build_normal_mode_matrix(p: usize) -> Vec<f64> {
    let pf = p as f64;
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut t = vec![0.0_f64; p * p];
    let inv_sqrt_p = 1.0 / pf.sqrt();
    let sqrt_2_over_p = (2.0 / pf).sqrt();
    for n in 0..p {
        t[n] = inv_sqrt_p;
    }
    let half = p / 2;
    for k in 1..=half {
        for n in 0..p {
            let angle = two_pi * (k as f64) * (n as f64) / pf;
            if p.is_multiple_of(2) && k == half {
                t[k * p + n] = inv_sqrt_p * if n % 2 == 0 { 1.0 } else { -1.0 };
            } else {
                t[k * p + n] = sqrt_2_over_p * angle.cos();
            }
        }
    }
    let sin_max = if p.is_multiple_of(2) { half - 1 } else { half };
    for k in 1..=sin_max {
        let row = p - k;
        for n in 0..p {
            let angle = two_pi * (k as f64) * (n as f64) / pf;
            t[row * p + n] = sqrt_2_over_p * angle.sin();
        }
    }
    t
}
/// Staging-coordinate transform for a single set of replica positions.
///
/// Returns staging vectors u\[k\]:
/// u\[0\] = r\[0\]
/// u\[k\] = r\[k\] - (k/(k+1)) * r\[k+1\] - (1/(k+1)) * r\[0\]  for k = 1..P-1
pub fn normal_mode_transform(replica_positions: &[[f64; 3]], n_beads: usize) -> Vec<[f64; 3]> {
    if n_beads == 0 {
        return vec![];
    }
    if n_beads == 1 {
        return replica_positions.to_vec();
    }
    let p = n_beads;
    let mut u = vec![[0.0f64; 3]; p];
    u[0] = replica_positions[0];
    for k in 1..p {
        let frac_next = k as f64 / (k as f64 + 1.0);
        let frac_first = 1.0 / (k as f64 + 1.0);
        let next = if k + 1 < p {
            replica_positions[k + 1]
        } else {
            replica_positions[0]
        };
        for d in 0..3 {
            u[k][d] = replica_positions[k][d]
                - frac_next * next[d]
                - frac_first * replica_positions[0][d];
        }
    }
    u
}
/// Staging fictitious mass factor for bead `k` of `n_beads`:
///
/// m_factor(k) = k * (n_beads - k + 1) / (n_beads - k + 2)
///
/// This is the relative factor for the fictitious staging mass (unitless).
pub fn staging_mass_factor(k: usize, n_beads: usize) -> f64 {
    if k == 0 || n_beads == 0 {
        return 1.0;
    }
    let kf = k as f64;
    let p = n_beads as f64;
    kf * (p - kf + 1.0) / (p - kf + 2.0)
}
/// Primitive kinetic energy estimator for PIMD.
///
/// KE_prim = (3*N*P/2)*k_B*T - spring_potential
///
/// where `spring_potential` is the sum of all inter-bead harmonic energies.
pub fn primitive_kinetic_estimator(
    t: f64,
    n_beads: usize,
    n_atoms: usize,
    spring_potential: f64,
) -> f64 {
    1.5 * n_atoms as f64 * n_beads as f64 * KB * t - spring_potential
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    fn approx_eq3(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
        (0..3).all(|d| (a[d] - b[d]).abs() < tol)
    }
    #[test]
    fn test_ring_polymer_centroid() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.positions[0] = [1.0, 0.0, 0.0];
        ring.positions[1] = [0.0, 1.0, 0.0];
        ring.positions[2] = [-1.0, 0.0, 0.0];
        ring.positions[3] = [0.0, -1.0, 0.0];
        let c = ring.centroid();
        assert!(approx_eq3(c, [0.0, 0.0, 0.0], 1e-15), "centroid = {:?}", c);
    }
    #[test]
    fn test_centroid_velocity() {
        let mut ring = RingPolymer::new(2, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.velocities[0] = [1.0, 2.0, 3.0];
        ring.velocities[1] = [3.0, 4.0, 5.0];
        let cv = ring.centroid_velocity();
        assert!(approx_eq3(cv, [2.0, 3.0, 4.0], 1e-12));
    }
    #[test]
    fn test_spring_forces_linear() {
        let mut ring = RingPolymer::new(2, 1.0e-27, [0.0, 0.0, 0.0]);
        let d = 1.0e-12_f64;
        ring.positions[0] = [d, 0.0, 0.0];
        ring.positions[1] = [-d, 0.0, 0.0];
        let temp = 300.0_f64;
        let forces = ring.spring_forces(temp);
        assert!(forces[0][0] < 0.0, "bead 0 force should be negative");
        assert!(forces[1][0] > 0.0, "bead 1 force should be positive");
        assert!(
            (forces[0][0] + forces[1][0]).abs() < 1e10,
            "forces should sum to ~zero"
        );
    }
    #[test]
    fn test_spring_potential_energy_positive() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.positions[0] = [1.0e-11, 0.0, 0.0];
        ring.positions[1] = [0.0, 1.0e-11, 0.0];
        ring.positions[2] = [-1.0e-11, 0.0, 0.0];
        ring.positions[3] = [0.0, -1.0e-11, 0.0];
        let pe = ring.spring_potential_energy(300.0);
        assert!(pe > 0.0, "spring PE should be positive for spread beads");
    }
    #[test]
    fn test_spring_potential_energy_zero_for_collapsed() {
        let ring = RingPolymer::new(4, 1.0e-27, [1.0e-10, 0.0, 0.0]);
        let pe = ring.spring_potential_energy(300.0);
        assert!(
            pe.abs() < 1e-50,
            "spring PE should be ~0 for collapsed ring"
        );
    }
    #[test]
    fn test_radius_of_gyration() {
        let mut ring_spread = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        ring_spread.positions[0] = [1.0, 0.0, 0.0];
        ring_spread.positions[1] = [0.0, 1.0, 0.0];
        ring_spread.positions[2] = [-1.0, 0.0, 0.0];
        ring_spread.positions[3] = [0.0, -1.0, 0.0];
        let ring_collapsed = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        assert!(ring_spread.radius_of_gyration() > ring_collapsed.radius_of_gyration());
        assert!(ring_collapsed.radius_of_gyration().abs() < 1e-15);
    }
    #[test]
    fn test_staging_roundtrip() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.positions[0] = [1.0, 2.0, 3.0];
        ring.positions[1] = [4.0, 5.0, 6.0];
        ring.positions[2] = [7.0, 8.0, 9.0];
        ring.positions[3] = [10.0, 11.0, 12.0];
        let staging = ring.to_staging();
        let recovered = RingPolymer::from_staging(&staging, ring.n_beads);
        for i in 0..ring.n_beads {
            assert!(
                approx_eq3(recovered[i], ring.positions[i], 1e-10),
                "bead {i}: {:?} != {:?}",
                recovered[i],
                ring.positions[i]
            );
        }
    }
    #[test]
    fn test_normal_mode_roundtrip() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.positions[0] = [1.0, 2.0, 3.0];
        ring.positions[1] = [4.0, 5.0, 6.0];
        ring.positions[2] = [7.0, 8.0, 9.0];
        ring.positions[3] = [10.0, 11.0, 12.0];
        let modes = ring.to_normal_modes();
        let recovered = RingPolymer::from_normal_modes(&modes, ring.n_beads);
        for i in 0..ring.n_beads {
            assert!(
                approx_eq3(recovered[i], ring.positions[i], 1e-10),
                "bead {i}: {:?} != {:?}",
                recovered[i],
                ring.positions[i]
            );
        }
    }
    #[test]
    fn test_normal_mode_centroid_is_centroid() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.positions[0] = [1.0, 0.0, 0.0];
        ring.positions[1] = [0.0, 1.0, 0.0];
        ring.positions[2] = [-1.0, 0.0, 0.0];
        ring.positions[3] = [0.0, -1.0, 0.0];
        let modes = ring.to_normal_modes();
        let centroid = ring.centroid();
        assert!(approx_eq3(modes[0], centroid, 1e-12));
    }
    #[test]
    fn test_pimd_step_conserves_topology() {
        let n_beads = 6_usize;
        let mut ring = RingPolymer::new(n_beads, 1.0e-27, [1.0e-10, 0.0, 0.0]);
        for (i, pos) in ring.positions.iter_mut().enumerate() {
            let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n_beads as f64);
            pos[0] += 1.0e-11 * angle.cos();
            pos[1] += 1.0e-11 * angle.sin();
        }
        PimdStep::velocity_verlet(&mut ring, [0.0; 3], 300.0, 1.0e-18);
        assert_eq!(ring.positions.len(), n_beads);
        assert_eq!(ring.velocities.len(), n_beads);
        for i in 0..n_beads {
            for d in 0..3 {
                assert!(ring.positions[i][d].is_finite());
                assert!(ring.velocities[i][d].is_finite());
            }
        }
    }
    #[test]
    fn test_primitive_kinetic_energy_classical_limit() {
        let ring = RingPolymer::new(1, 1.0e-27, [1.0e-10, 0.0, 0.0]);
        let temp = 300.0;
        let ke = ring.primitive_kinetic_energy(temp);
        let expected = 1.5 * KB * temp;
        assert!(
            (ke - expected).abs() < 1e-30,
            "single-bead KE = {ke}, expected {expected}"
        );
    }
    #[test]
    fn test_virial_kinetic_energy() {
        let ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        let ext_forces = vec![[0.0; 3]; 4];
        let ke = ring.virial_kinetic_energy(300.0, &ext_forces);
        let expected = 1.5 * KB * 300.0;
        assert!(
            (ke - expected).abs() < 1e-30,
            "virial KE = {ke}, expected {expected}"
        );
    }
    #[test]
    fn test_potential_energy_estimator() {
        let ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        let potentials = vec![1.0, 2.0, 3.0, 4.0];
        let v = ring.potential_energy_estimator(&potentials);
        assert!((v - 2.5).abs() < 1e-12, "V_est = {v}, expected 2.5");
    }
    #[test]
    fn test_classical_kinetic_energy() {
        let mut ring = RingPolymer::new(2, 1.0e-27, [0.0, 0.0, 0.0]);
        ring.velocities[0] = [1.0, 0.0, 0.0];
        ring.velocities[1] = [0.0, 1.0, 0.0];
        let ke = ring.classical_kinetic_energy();
        let expected = 0.5 * 1.0e-27 * 1.0 + 0.5 * 1.0e-27 * 1.0;
        assert!(
            (ke - expected).abs() < 1e-40,
            "KE = {ke}, expected {expected}"
        );
    }
    #[test]
    fn test_rpmd_hamiltonian_finite() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [1.0e-10, 0.0, 0.0]);
        for (i, pos) in ring.positions.iter_mut().enumerate() {
            pos[0] += (i as f64) * 1e-12;
        }
        ring.velocities[0] = [100.0, 0.0, 0.0];
        let potentials = vec![0.0; 4];
        let h = RpmdStep::hamiltonian(&ring, 300.0, &potentials);
        assert!(h.is_finite(), "Hamiltonian should be finite");
        assert!(h >= 0.0, "Hamiltonian should be non-negative");
    }
    #[test]
    fn test_pile_centroid_thermostat() {
        let pile = PileThermostat::new(1e12, 300.0);
        let mut cv = [1000.0, 0.0, 0.0];
        let noise = [0.0, 0.0, 0.0];
        pile.apply_centroid(&mut cv, 1.0e-27, 1e-15, noise);
        assert!(cv[0].abs() < 1000.0, "velocity should be damped");
    }
    #[test]
    fn test_pile_optimal_gamma() {
        let pile = PileThermostat::new(1e12, 300.0);
        let gamma = pile.optimal_gamma_internal(4);
        let expected = 4.0 * KB * 300.0 / HBAR;
        assert!((gamma - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_bead_propagator_frequencies() {
        let freqs = BeadPropagator::normal_mode_frequencies(4, 1.0e-27, 300.0);
        assert_eq!(freqs.len(), 4);
        assert!(freqs[0].abs() < 1e-10, "mode 0 freq should be ~0");
        assert!(freqs[1] > 0.0);
    }
    #[test]
    fn test_bead_propagator_free_particle() {
        let mut q = [1.0, 0.0, 0.0];
        let mut v = [2.0, 0.0, 0.0];
        let dt = 0.5;
        BeadPropagator::propagate_mode(&mut q, &mut v, 0.0, dt);
        assert!((q[0] - 2.0).abs() < 1e-12, "q = {}, expected 2.0", q[0]);
        assert!((v[0] - 2.0).abs() < 1e-12, "v unchanged for free particle");
    }
    #[test]
    fn test_bead_propagator_harmonic_energy_conservation() {
        let omega = 1.0;
        let mut q = [1.0, 0.0, 0.0];
        let mut v = [0.0, 0.0, 0.0];
        let dt = 0.01;
        let mass = 1.0;
        let e0 = 0.5 * mass * v[0] * v[0] + 0.5 * mass * omega * omega * q[0] * q[0];
        for _ in 0..1000 {
            BeadPropagator::propagate_mode(&mut q, &mut v, omega, dt);
        }
        let e1 = 0.5 * mass * v[0] * v[0] + 0.5 * mass * omega * omega * q[0] * q[0];
        assert!(
            (e1 - e0).abs() / e0 < 1e-10,
            "energy not conserved: {e0} -> {e1}"
        );
    }
    #[test]
    fn test_rpmd_multi_step() {
        let n = 4;
        let mut ring = RingPolymer::new(n, 1.0e-27, [1.0e-10, 0.0, 0.0]);
        let forces = vec![[0.0; 3]; n];
        RpmdStep::multi_step(&mut ring, &forces, 300.0, 1e-18, 10);
        for i in 0..n {
            for d in 0..3 {
                assert!(ring.positions[i][d].is_finite());
            }
        }
    }
    #[test]
    fn test_centroid_md_step_moves_centroid() {
        let n = 4;
        let mut ring = RingPolymer::new(n, 1.0e-27, [0.0, 0.0, 0.0]);
        let f_ext = 1.0e-12;
        let per_bead_forces: Vec<[f64; 3]> = vec![[f_ext, 0.0, 0.0]; n];
        let c0 = ring.centroid();
        CentroidMD::step(&mut ring, &per_bead_forces, 1e-15);
        let c1 = ring.centroid();
        assert!(
            c1[0].abs() > c0[0].abs() || (c1[0] - c0[0]).abs() > 0.0,
            "centroid should move after CMD step"
        );
    }
    #[test]
    fn test_centroid_md_msd_zero_initial() {
        let ring = RingPolymer::new(4, 1.0e-27, [1.0e-10, 0.0, 0.0]);
        let init_c = ring.centroid();
        let msd = CentroidMD::centroid_msd(&ring, init_c);
        assert!(msd.abs() < 1e-40, "MSD from initial position should be 0");
    }
    #[test]
    fn test_centroid_md_msd_after_move() {
        let mut ring = RingPolymer::new(4, 1.0e-27, [0.0, 0.0, 0.0]);
        let init_c = ring.centroid();
        for pos in ring.positions.iter_mut() {
            pos[0] += 1.0e-10;
        }
        let msd = CentroidMD::centroid_msd(&ring, init_c);
        assert!(msd > 0.0, "MSD should be positive after moving beads");
    }
    #[test]
    fn test_feynman_hibbs_prefactor_positive() {
        let fh = FeynmanHibbsCorrection::new(1.674e-27, 300.0);
        let pf = fh.prefactor();
        assert!(
            pf > 0.0,
            "Feynman-Hibbs prefactor should be positive, got {pf}"
        );
    }
    #[test]
    fn test_feynman_hibbs_effective_potential_with_zero_curvature() {
        let fh = FeynmanHibbsCorrection::new(1.674e-27, 300.0);
        let v = 1.0e-21;
        let v_eff = fh.effective_potential(v, 0.0);
        assert!((v_eff - v).abs() < 1e-40, "with zero curvature V_eff = V");
    }
    #[test]
    fn test_feynman_hibbs_effective_potential_correction_sign() {
        let fh = FeynmanHibbsCorrection::new(1.674e-27, 300.0);
        let v_eff = fh.effective_potential(0.0, 1e40);
        assert!(v_eff > 0.0, "positive curvature should increase V_eff");
    }
    #[test]
    fn test_thermal_de_broglie_wavelength_positive() {
        let fh = FeynmanHibbsCorrection::new(1.674e-27, 300.0);
        let lam = fh.thermal_de_broglie_wavelength();
        assert!(
            lam > 0.0 && lam.is_finite(),
            "thermal de Broglie wavelength = {lam}"
        );
    }
    #[test]
    fn test_wkb_tunneling_probability_bounded() {
        let fh = FeynmanHibbsCorrection::new(9.109e-31, 300.0);
        let p = fh.wkb_tunneling_probability(1.0e-19, 1.0e-10);
        assert!(
            (0.0..=1.0).contains(&p),
            "WKB probability out of [0,1]: {p}"
        );
    }
    #[test]
    fn test_wkb_tunneling_zero_barrier() {
        let fh = FeynmanHibbsCorrection::new(9.109e-31, 300.0);
        let p = fh.wkb_tunneling_probability(0.0, 1.0e-10);
        assert!((p - 1.0).abs() < 1e-10, "zero barrier → T = 1, got {p}");
    }
    #[test]
    fn test_staging_fictitious_masses_length() {
        let masses = StagingTransform::fictitious_masses(4, 1.0e-27);
        assert_eq!(masses.len(), 4);
    }
    #[test]
    fn test_staging_fictitious_masses_values() {
        let m = 1.0e-27;
        let masses = StagingTransform::fictitious_masses(4, m);
        assert!((masses[0] - m).abs() < 1e-50, "m_0 should be m");
        assert!((masses[1] - 2.0 * m).abs() < 1e-50, "m_1 = 2m");
        assert!((masses[2] - 1.5 * m).abs() < 1e-50, "m_2 = 1.5m");
    }
    #[test]
    fn test_staging_frequencies_centroid_zero() {
        let freqs = StagingTransform::staging_frequencies(4, 300.0);
        assert!(freqs[0].abs() < 1e-50, "centroid frequency should be 0");
    }
    #[test]
    fn test_staging_frequencies_positive_for_internal() {
        let freqs = StagingTransform::staging_frequencies(4, 300.0);
        for k in 1..4 {
            assert!(freqs[k] > 0.0, "mode {k} frequency should be positive");
        }
    }
    #[test]
    fn test_staging_proper_first_bead_preserved() {
        let positions = vec![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ];
        let u = StagingTransform::to_staging_proper(&positions, 4);
        for d in 0..3 {
            assert!(
                (u[0][d] - positions[0][d]).abs() < 1e-12,
                "u[0][{d}] = {}, expected {}",
                u[0][d],
                positions[0][d]
            );
        }
    }
    #[test]
    fn test_pimd_thermostat_gamma_lengths() {
        let therm = PimdThermostat::new_optimal(4, 1.0e-27, 300.0);
        assert_eq!(therm.gamma.len(), 4);
    }
    #[test]
    fn test_pimd_thermostat_centroid_gamma_zero() {
        let therm = PimdThermostat::new_optimal(4, 1.0e-27, 300.0);
        assert!(therm.gamma[0].abs() < 1e-50, "centroid gamma should be 0");
    }
    #[test]
    fn test_pimd_thermostat_internal_gamma_positive() {
        let therm = PimdThermostat::new_optimal(4, 1.0e-27, 300.0);
        for k in 1..4 {
            assert!(therm.gamma[k] > 0.0, "mode {k} gamma should be positive");
        }
    }
    #[test]
    fn test_pimd_thermostat_apply_mode_damps_velocity() {
        let therm = PimdThermostat::new_optimal(4, 1.0e-27, 300.0);
        let mut vel = [1000.0, 0.0, 0.0];
        let noise = [0.0, 0.0, 0.0];
        therm.apply_mode(&mut vel, 2, 1e-15, noise);
        assert!(vel[0].abs() < 1000.0, "velocity should be damped");
    }
    #[test]
    fn test_pimd_atom_centroid_equal_replicas() {
        let n_beads = 4;
        let mut atom = PimdAtom::new(n_beads, 1.0e-27);
        for r in atom.replicas.iter_mut() {
            *r = [2.0, 3.0, 4.0];
        }
        let c = atom.centroid();
        assert!((c[0] - 2.0).abs() < 1e-12);
        assert!((c[1] - 3.0).abs() < 1e-12);
        assert!((c[2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_pimd_atom_spring_force_zero_when_aligned() {
        let n_beads = 4;
        let mut atom = PimdAtom::new(n_beads, 1.0e-27);
        for r in atom.replicas.iter_mut() {
            *r = [1.0, 0.0, 0.0];
        }
        let fs = atom.spring_forces(300.0);
        for f in &fs {
            for d in 0..3 {
                assert!(
                    f[d].abs() < 1e-30,
                    "spring force should be zero for aligned replicas"
                );
            }
        }
    }
    #[test]
    fn test_pimd_state_quantum_ke_positive() {
        let n_beads = 4;
        let n_atoms = 2;
        let temperature = 300.0;
        let mut state = PimdState {
            atoms: (0..n_atoms)
                .map(|_| PimdAtom::new(n_beads, 1.0e-27))
                .collect(),
            temperature,
            n_beads,
        };
        for atom in state.atoms.iter_mut() {
            for r in atom.replicas.iter_mut() {
                *r = [0.0; 3];
            }
        }
        let ke = state.quantum_kinetic_energy();
        assert!(ke > 0.0, "quantum KE should be positive, got {ke}");
    }
    #[test]
    fn test_normal_mode_transform_lossless() {
        let positions: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let u = normal_mode_transform(&positions, 4);
        assert!((u[0][0] - 1.0).abs() < 1e-12, "u[0] should equal r[0]");
        for k in 1..4 {
            for d in 0..3 {
                assert!(
                    u[k][d].abs() < 1e-12,
                    "u[{k}][{d}] should be ~0 for equal replicas, got {}",
                    u[k][d]
                );
            }
        }
    }
}
/// Full ring-polymer normal-mode transformation matrix.
///
/// Returns the P×P orthogonal matrix C where the k-th normal mode coordinate
/// is ũ_k = Σ_j C_{kj} r_j.  The zeroth mode is the centroid.
pub fn normal_mode_matrix(p: usize) -> Vec<Vec<f64>> {
    use std::f64::consts::PI;
    let mut c = vec![vec![0.0f64; p]; p];
    for j in 0..p {
        c[0][j] = 1.0 / (p as f64).sqrt();
        for k in 1..=p / 2 {
            let arg = 2.0 * PI * k as f64 * j as f64 / p as f64;
            if k < p - k {
                c[k][j] = (2.0 / p as f64).sqrt() * arg.cos();
                c[p - k][j] = (2.0 / p as f64).sqrt() * arg.sin();
            } else {
                c[k][j] = (1.0 / p as f64).sqrt() * if j % 2 == 0 { 1.0 } else { -1.0 };
            }
        }
    }
    c
}
/// Transform bead positions to normal-mode coordinates using the full matrix.
pub fn beads_to_normal_modes(beads: &[[f64; 3]], c: &[Vec<f64>]) -> Vec<[f64; 3]> {
    let p = beads.len();
    let mut modes = vec![[0.0f64; 3]; p];
    for k in 0..p {
        for j in 0..p {
            for d in 0..3 {
                modes[k][d] += c[k][j] * beads[j][d];
            }
        }
    }
    modes
}
/// Transform normal-mode coordinates back to bead positions (inverse = transpose for orthogonal C).
pub fn normal_modes_to_beads(modes: &[[f64; 3]], c: &[Vec<f64>]) -> Vec<[f64; 3]> {
    let p = modes.len();
    let mut beads = vec![[0.0f64; 3]; p];
    for j in 0..p {
        for k in 0..p {
            for d in 0..3 {
                beads[j][d] += c[k][j] * modes[k][d];
            }
        }
    }
    beads
}
/// Normal-mode spring frequencies ω_k for a ring polymer of P beads at temperature T.
///
/// ω_k = 2 ω_P sin(k π / P),  ω_P = P k_B T / ħ.
pub fn normal_mode_frequencies(p: usize, temperature: f64) -> Vec<f64> {
    let omega_p = p as f64 * KB * temperature / HBAR;
    let mut freqs = vec![0.0f64; p];
    for k in 0..p {
        freqs[k] = 2.0 * omega_p * (k as f64 * PI / p as f64).sin().abs();
    }
    freqs
}
/// Quantum partition function Z(β) for a harmonic oscillator.
///
/// Z = 1 / (2 sinh(βħω/2)) where β = 1/(k_B T).
pub fn harmonic_partition_function(omega: f64, temperature: f64) -> f64 {
    let beta = 1.0 / (KB * temperature);
    let x = 0.5 * beta * HBAR * omega;
    1.0 / (2.0 * x.sinh())
}
/// Zero-point energy of a harmonic oscillator (J).
///
/// E_ZPE = ħ ω / 2.
pub fn zero_point_energy(omega: f64) -> f64 {
    0.5 * HBAR * omega
}
/// Thermal de Broglie wavelength (m).
///
/// Λ = ħ sqrt(2π β / m).
pub fn thermal_de_broglie(mass: f64, temperature: f64) -> f64 {
    let beta = 1.0 / (KB * temperature);
    HBAR * (2.0 * std::f64::consts::PI * beta / mass).sqrt()
}
/// Feynman–Hibbs effective potential correction.
///
/// To lowest order in ħ², the quantum effective potential is:
/// V_eff(x) ≈ V(x) + (ħ²β/24m) V''(x)
///
/// `v_second_deriv` is ∂²V/∂x² evaluated at x.
pub fn feynman_hibbs_correction(v: f64, v_second_deriv: f64, mass: f64, temperature: f64) -> f64 {
    let beta = 1.0 / (KB * temperature);
    let correction = HBAR * HBAR * beta / (24.0 * mass) * v_second_deriv;
    v + correction
}
/// WKB tunnelling probability through a parabolic barrier.
///
/// P_tunnel = exp(-2 S / ħ) where S = π m ω_b Δx / 2
/// for a parabolic barrier of height V_b, width 2a, with ω_b = sqrt(2 V_b / (m a²)).
///
/// Returns the WKB tunnelling exponent 2S/ħ (dimensionless).
pub fn wkb_tunnel_exponent(v_barrier: f64, half_width: f64, mass: f64) -> f64 {
    let omega_b = (2.0 * v_barrier / (mass * half_width * half_width)).sqrt();
    let s = std::f64::consts::PI * mass * omega_b * half_width * half_width / 2.0;
    2.0 * s / HBAR
}
/// Tunnel splitting estimate from instanton theory.
///
/// Δ ≈ ħ ω₀ / π * exp(-S_inst / ħ) where S_inst is the instanton action.
///
/// For a double-well of height V_b, frequency ω₀ at the minima, the instanton
/// action is approximately S_inst = π m ω₀ a² / 2 (same as WKB for parabolic).
pub fn tunnel_splitting(omega0: f64, s_inst_over_hbar: f64) -> f64 {
    HBAR * omega0 / std::f64::consts::PI * (-s_inst_over_hbar).exp()
}
/// Virial kinetic energy estimator for path integral simulations.
///
/// KE_vir = (3N/2) k_B T + (1/2P) Σ_{i,s} (r_{is} - r̄_i) · F_{is}
///
/// where r̄_i is the centroid of atom i and F_{is} is the physical force on
/// bead s of atom i.
///
/// This estimator has lower variance than the primitive estimator.
pub fn virial_kinetic_estimator(
    beads: &[Vec<[f64; 3]>],
    forces: &[Vec<[f64; 3]>],
    temperature: f64,
    n_beads: usize,
) -> f64 {
    let n_atoms = beads.len();
    let classical_part = 1.5 * n_atoms as f64 * KB * temperature;
    let mut quantum_correction = 0.0f64;
    for i in 0..n_atoms {
        let mut centroid = [0.0f64; 3];
        for s in 0..n_beads {
            for d in 0..3 {
                centroid[d] += beads[i][s][d];
            }
        }
        for d in 0..3 {
            centroid[d] /= n_beads as f64;
        }
        for s in 0..n_beads {
            for d in 0..3 {
                quantum_correction += (beads[i][s][d] - centroid[d]) * forces[i][s][d];
            }
        }
    }
    classical_part + 0.5 / n_beads as f64 * quantum_correction
}
#[cfg(test)]
mod tests_pimd_ext {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_normal_mode_matrix_orthogonal() {
        let p = 4;
        let c = normal_mode_matrix(p);
        for i in 0..p {
            for j in 0..p {
                let dot: f64 = (0..p).map(|k| c[i][k] * c[j][k]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-12,
                    "C C^T [{i}][{j}] = {dot}, expected {expected}"
                );
            }
        }
    }
    #[test]
    fn test_beads_to_normal_modes_roundtrip() {
        let p = 4;
        let c = normal_mode_matrix(p);
        let original: Vec<[f64; 3]> = vec![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ];
        let modes = beads_to_normal_modes(&original, &c);
        let recovered = normal_modes_to_beads(&modes, &c);
        for i in 0..p {
            for d in 0..3 {
                assert!(
                    (recovered[i][d] - original[i][d]).abs() < 1e-10,
                    "roundtrip failed at bead {i} dim {d}"
                );
            }
        }
    }
    #[test]
    fn test_normal_mode_frequencies_centroid_zero() {
        let freqs = normal_mode_frequencies(4, 300.0);
        assert!(freqs[0].abs() < 1e-30, "centroid frequency should be zero");
    }
    #[test]
    fn test_normal_mode_frequencies_positive() {
        let freqs = normal_mode_frequencies(4, 300.0);
        for k in 1..4 {
            assert!(freqs[k] > 0.0, "mode {k} frequency should be positive");
        }
    }
    #[test]
    fn test_harmonic_partition_function_high_t() {
        let omega = 1e12_f64;
        let t = 10000.0;
        let z = harmonic_partition_function(omega, t);
        let z_classical = KB * t / (HBAR * omega);
        let rel_err = (z - z_classical).abs() / z_classical;
        assert!(
            rel_err < 0.01,
            "high-T partition function rel error = {rel_err}"
        );
    }
    #[test]
    fn test_zero_point_energy_positive() {
        let zpe = zero_point_energy(1e12);
        assert!(zpe > 0.0, "zero-point energy should be positive");
    }
    #[test]
    fn test_thermal_de_broglie_decreases_with_temperature() {
        let m = 1.67e-27;
        let lambda_cold = thermal_de_broglie(m, 10.0);
        let lambda_hot = thermal_de_broglie(m, 300.0);
        assert!(
            lambda_cold > lambda_hot,
            "de Broglie wavelength should decrease with T"
        );
    }
    #[test]
    fn test_feynman_hibbs_correction_positive_for_positive_curvature() {
        let v = 0.0;
        let v_pp = 1.0e40;
        let corr = feynman_hibbs_correction(v, v_pp, 1e-27, 300.0);
        assert!(
            corr > v,
            "Feynman-Hibbs correction should be positive for positive V''"
        );
    }
    #[test]
    fn test_wkb_tunnel_exponent_positive() {
        let exp = wkb_tunnel_exponent(1.0e-20, 1.0e-10, 1.67e-27);
        assert!(exp > 0.0, "WKB exponent should be positive");
    }
    #[test]
    fn test_tunnel_splitting_decreases_with_exponent() {
        let delta1 = tunnel_splitting(1e12, 10.0);
        let delta2 = tunnel_splitting(1e12, 20.0);
        assert!(
            delta1 > delta2,
            "tunnel splitting should decrease with larger exponent"
        );
    }
    #[test]
    fn test_rpc_contraction_preserves_centroid() {
        let p_full = 8;
        let p_cont = 4;
        let rpc = RingPolymerContraction::new(p_full, p_cont);
        let beads: Vec<[f64; 3]> = (0..p_full).map(|_| [1.0, 2.0, 3.0]).collect();
        let contracted = rpc.contract(&beads, 300.0);
        let mut c = [0.0f64; 3];
        for b in &contracted {
            for d in 0..3 {
                c[d] += b[d];
            }
        }
        for d in 0..3 {
            c[d] /= p_cont as f64;
        }
        assert!(
            (c[0] - 1.0).abs() < 1e-10,
            "centroid x should be 1.0, got {}",
            c[0]
        );
    }
    #[test]
    fn test_virial_estimator_pure_classical() {
        let p = 4;
        let n = 2;
        let temp = 300.0;
        let beads: Vec<Vec<[f64; 3]>> = (0..n).map(|_| vec![[0.0, 0.0, 0.0]; p]).collect();
        let forces: Vec<Vec<[f64; 3]>> = (0..n).map(|_| vec![[1.0, 2.0, 3.0]; p]).collect();
        let ke = virial_kinetic_estimator(&beads, &forces, temp, p);
        let ke_classical = 1.5 * n as f64 * KB * temp;
        assert!(
            (ke - ke_classical).abs() < 1e-40,
            "virial KE should equal classical part when beads coincide"
        );
    }
    #[test]
    fn test_pile_l_thermostat_centroid_friction() {
        let th = PileLThermostat::new(4, 300.0, 1.0);
        assert!(
            (th.gamma[0] - 1.0).abs() < 1e-15,
            "centroid gamma should equal gamma0"
        );
    }
    #[test]
    fn test_pile_l_thermostat_higher_modes_critical_damping() {
        let th = PileLThermostat::new(4, 300.0, 1.0);
        for k in 1..4 {
            let expected = 2.0 * th.omega[k];
            assert!(
                (th.gamma[k] - expected).abs() < 1e-10,
                "mode {k} gamma should be 2*omega, got {} vs {expected}",
                th.gamma[k]
            );
        }
    }
    #[test]
    fn test_trpmd_integrator_builds() {
        let integrator = TrpmdIntegrator::new(4, 300.0, 1.0);
        assert_eq!(integrator.n_beads, 4);
    }
    #[test]
    fn test_normal_mode_propagator_centroid_free_motion() {
        let prop = NormalModePropagator::new(4, 300.0);
        let pos: Vec<[f64; 3]> = vec![[0.0; 3]; 4];
        let vel: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; 4];
        let mass = 1.67e-27;
        let dt = 1e-15;
        let (new_pos, _new_vel) = prop.propagate(&pos, &vel, mass, dt);
        let centroid_x: f64 = new_pos.iter().map(|p| p[0]).sum::<f64>() / 4.0;
        assert!(
            (centroid_x - dt).abs() < 1e-17,
            "centroid should move by v*dt = {dt}, got {centroid_x}"
        );
    }
}
