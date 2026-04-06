//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Reduced Planck constant (J·s).
pub(super) const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant (J/K).
pub(super) const KB: f64 = 1.380_649e-23;
/// Gas constant (J mol⁻¹ K⁻¹).
pub(super) const R_GAS: f64 = 8.314_462_618;
/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Scale a 3-vector.
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Add two 3-vectors.
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Norm of a 3-vector.
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Linear interpolation between two 3-vectors.
pub(super) fn interpolate3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    add3(scale3(a, 1.0 - t), scale3(b, t))
}
/// Assemble a simplified 3×3 mass-weighted dynamical matrix for a mono-atomic
/// crystal with isotropic harmonic spring constant `k`.
///
/// Returns a flat 3×3 matrix (row-major).
pub fn dynamical_matrix(k_spring: f64, mass: f64, coordination: usize) -> [f64; 9] {
    let d = k_spring * coordination as f64 / mass;
    [d, 0.0, 0.0, 0.0, d, 0.0, 0.0, 0.0, d]
}
/// Extract phonon frequencies (in rad/s) from a 3×3 dynamical matrix via the
/// eigenvalues of the diagonal approximation.
pub fn phonon_frequencies(dyn_mat: &[f64; 9]) -> [f64; 3] {
    [
        dyn_mat[0].max(0.0).sqrt(),
        dyn_mat[4].max(0.0).sqrt(),
        dyn_mat[8].max(0.0).sqrt(),
    ]
}
/// Compute the phonon density of states (DOS) as a histogram over frequency.
///
/// Returns a `Vec`f64` of length `n_bins` with counts normalised to total 1.
pub fn density_of_states_phonon(frequencies: &[f64], n_bins: usize, omega_max: f64) -> Vec<f64> {
    let mut hist = vec![0usize; n_bins];
    let dw = omega_max / n_bins as f64;
    for &w in frequencies {
        if w >= 0.0 && w < omega_max {
            let bin = (w / dw).floor() as usize;
            let bin = bin.min(n_bins - 1);
            hist[bin] += 1;
        }
    }
    let total = frequencies.len() as f64;
    if total > 0.0 {
        hist.iter().map(|&c| c as f64 / total).collect()
    } else {
        vec![0.0; n_bins]
    }
}
/// Estimate the Debye temperature from the maximum phonon frequency.
///
/// `θ_D = ħ ω_max / k_B`
pub fn debye_temperature(omega_max: f64, hbar: f64, k_b: f64) -> f64 {
    hbar * omega_max / k_b
}
/// Einstein model molar heat capacity.
///
/// `C_V = 3R (θ_E/T)² exp(θ_E/T) / (exp(θ_E/T) − 1)²`
pub fn heat_capacity_einstein(theta_e: f64, temp: f64, r_gas: f64) -> f64 {
    if temp <= 0.0 {
        return 0.0;
    }
    let x = theta_e / temp;
    if x > 700.0 {
        return 0.0;
    }
    let ex = x.exp();
    let denom = (ex - 1.0).powi(2);
    if denom < 1e-300 {
        return 0.0;
    }
    3.0 * r_gas * x * x * ex / denom
}
/// Debye frequency from number density and spring constant (1D estimate).
///
/// ω_D = π / a * sqrt(k / m)  (zone-boundary frequency for 1D chain)
pub fn debye_frequency(n: f64, mass: f64, spring_const: f64) -> f64 {
    if mass <= 0.0 || n <= 0.0 {
        return 0.0;
    }
    let a = 1.0 / n;
    (PI / a) * (spring_const / mass).sqrt() / 2.0
}
/// Einstein frequency for an atom in a harmonic potential well.
///
/// ω_E = sqrt(k / m)
pub fn einstein_frequency(spring_const: f64, mass: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    (spring_const / mass).sqrt()
}
/// 1D phonon density of states at frequency ω.
///
/// g(ω) = 2 / (π √(ω_max² − ω²))
///
/// Returns `0` outside [0, ω_max).
pub fn phonon_density_of_states_1d(omega: f64, omega_max: f64) -> f64 {
    if omega < 0.0 || omega >= omega_max || omega_max <= 0.0 {
        return 0.0;
    }
    let d = omega_max * omega_max - omega * omega;
    if d < f64::EPSILON {
        return 0.0;
    }
    2.0 / (PI * d.sqrt())
}
/// Bose-Einstein occupation number.
///
/// n(ω) = 1 / (exp(ℏω / kT) − 1)
///
/// Returns `0` at T = 0 and `f64::INFINITY` at ω = 0.
pub fn bose_einstein_occupation(omega: f64, temperature: f64, hbar: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    if omega <= 0.0 {
        return f64::INFINITY;
    }
    let x = hbar * omega / (KB * temperature);
    if x > 700.0 {
        return 0.0;
    }
    1.0 / (x.exp() - 1.0)
}
/// Phonon contribution to heat capacity summed over a set of modes.
///
/// C = Σ_i k_B x_i² exp(x_i) / (exp(x_i) − 1)²,  x_i = ℏω_i / kT
pub fn phonon_heat_capacity(omegas: &[f64], temperature: f64, hbar: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    omegas
        .iter()
        .map(|&w| {
            if w <= 0.0 {
                return 0.0;
            }
            let x = hbar * w / (KB * temperature);
            if x > 700.0 {
                return 0.0;
            }
            let ex = x.exp();
            let denom = (ex - 1.0) * (ex - 1.0);
            if denom < 1e-300 {
                return 0.0;
            }
            KB * x * x * ex / denom
        })
        .sum()
}
/// Grüneisen parameter γ = V β K_T / C_V.
///
/// # Arguments
/// * `volume_expansion` – Volumetric thermal expansion coefficient β (1/K).
/// * `bulk_modulus`     – Isothermal bulk modulus K_T (Pa or GPa).
/// * `heat_capacity`    – Constant-volume heat capacity C_V (J/K).
/// * `volume`           – Molar or unit-cell volume (m³).
pub fn gruneisen_parameter(
    volume_expansion: f64,
    bulk_modulus: f64,
    heat_capacity: f64,
    volume: f64,
) -> f64 {
    if heat_capacity.abs() < f64::EPSILON {
        return 0.0;
    }
    volume * volume_expansion * bulk_modulus / heat_capacity
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_dot3_orthogonal() {
        assert_eq!(dot3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 0.0);
    }
    #[test]
    fn test_dot3_parallel() {
        assert!((dot3([2.0, 0.0, 0.0], [2.0, 0.0, 0.0]) - 4.0).abs() < 1e-14);
    }
    #[test]
    fn test_cross3_standard() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-14);
        assert!((c[1]).abs() < 1e-14);
        assert!((c[2] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_cross3_anti_commute() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c1 = cross3(a, b);
        let c2 = cross3(b, a);
        for i in 0..3 {
            assert!((c1[i] + c2[i]).abs() < 1e-13);
        }
    }
    #[test]
    fn test_cubic_volume() {
        let lv = LatticeVectors::cubic(2.0);
        assert!((lv.volume() - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_cubic_volume_unit() {
        let lv = LatticeVectors::cubic(1.0);
        assert!((lv.volume() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_reciprocal_vectors_cubic() {
        let a = 2.0_f64;
        let lv = LatticeVectors::cubic(a);
        let (b1, b2, b3) = lv.reciprocal_vectors();
        assert!((b1[0] - 2.0 * PI / a).abs() < 1e-12);
        assert!(b1[1].abs() < 1e-12);
        assert!(b1[2].abs() < 1e-12);
        assert!((b2[1] - 2.0 * PI / a).abs() < 1e-12);
        assert!((b3[2] - 2.0 * PI / a).abs() < 1e-12);
    }
    #[test]
    fn test_reciprocal_orthogonality_cubic() {
        let lv = LatticeVectors::cubic(3.0);
        let (b1, b2, b3) = lv.reciprocal_vectors();
        assert!((dot3(lv.a1, b1) - 2.0 * PI).abs() < 1e-10);
        assert!(dot3(lv.a1, b2).abs() < 1e-10);
        assert!(dot3(lv.a1, b3).abs() < 1e-10);
        assert!((dot3(lv.a2, b2) - 2.0 * PI).abs() < 1e-10);
        assert!((dot3(lv.a3, b3) - 2.0 * PI).abs() < 1e-10);
    }
    #[test]
    fn test_lattice_vectors_new() {
        let lv = LatticeVectors::new([1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]);
        assert!((lv.volume() - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_fcc_volume() {
        let a = 4.0_f64;
        let lv = LatticeVectors::fcc(a);
        let expected = a * a * a / 4.0;
        assert!(
            (lv.volume().abs() - expected).abs() < 1e-10,
            "vol={}",
            lv.volume()
        );
    }
    #[test]
    fn test_bcc_volume() {
        let a = 3.0_f64;
        let lv = LatticeVectors::bcc(a);
        let expected = a * a * a / 2.0;
        assert!((lv.volume().abs() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_dynamical_matrix_diagonal() {
        let dm = dynamical_matrix(1.0, 1.0, 6);
        assert_eq!(dm[1], 0.0);
        assert_eq!(dm[2], 0.0);
        assert_eq!(dm[3], 0.0);
        assert_eq!(dm[5], 0.0);
    }
    #[test]
    fn test_dynamical_matrix_diagonal_value() {
        let dm = dynamical_matrix(2.0, 4.0, 3);
        assert!((dm[0] - 1.5).abs() < 1e-12);
        assert!((dm[4] - 1.5).abs() < 1e-12);
        assert!((dm[8] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_dynmat_new_cubic_frequencies_positive() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(1.0, 1.0, 6, lv);
        let freqs = dm.phonon_frequencies_gamma();
        for &w in &freqs {
            assert!(w >= 0.0);
        }
    }
    #[test]
    fn test_dynmat_acoustic_at_gamma() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(1.0, 1.0, 6, lv);
        assert!(dm.acoustic_at_gamma());
    }
    #[test]
    fn test_dynmat_max_frequency_positive() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(4.0, 1.0, 6, lv);
        assert!(dm.max_frequency() > 0.0);
    }
    #[test]
    fn test_dynmat_at_q_zone_boundary_nonzero() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(1.0, 1.0, 6, lv);
        let freqs = dm.at_q([PI, 0.0, 0.0]);
        assert!(freqs[0] > 0.0);
    }
    #[test]
    fn test_dispersion_n_kpoints() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(1.0, 1.0, 6, lv);
        let disp = PhononDispersion::compute_sc(dm, 10);
        assert_eq!(disp.n_kpoints(), 31);
    }
    #[test]
    fn test_dispersion_max_frequency_positive() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(4.0, 1.0, 6, lv);
        let disp = PhononDispersion::compute_sc(dm, 20);
        assert!(disp.max_frequency() > 0.0);
    }
    #[test]
    fn test_dispersion_acoustic_at_gamma() {
        let lv = LatticeVectors::cubic(1.0);
        let dm = DynamicalMatrix::new_cubic(1.0, 1.0, 6, lv);
        let disp = PhononDispersion::compute_sc(dm, 10);
        assert!(disp.acoustic_modes_at_gamma(1e-6));
    }
    #[test]
    fn test_dispersion_all_frequencies_nonnegative() {
        let lv = LatticeVectors::cubic(2.0);
        let dm = DynamicalMatrix::new_cubic(1.0, 1.0, 6, lv);
        let disp = PhononDispersion::compute_sc(dm, 15);
        for &w in disp.all_frequencies().iter() {
            assert!(w >= 0.0);
        }
    }
    #[test]
    fn test_phonon_dos_norm_approx_one() {
        let freqs: Vec<f64> = (1..=100).map(|i| i as f64 * 0.1).collect();
        let dos = PhononDOS::from_frequencies(&freqs, 20);
        let n = dos.norm();
        assert!(n > 0.5 && n < 1.5, "norm = {n}");
    }
    #[test]
    fn test_phonon_dos_empty_input() {
        let dos = PhononDOS::from_frequencies(&[], 10);
        assert!(dos.dos.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_phonon_dos_peak_positive() {
        let freqs: Vec<f64> = (1..=50).map(|i| i as f64).collect();
        let dos = PhononDOS::from_frequencies(&freqs, 10);
        assert!(dos.peak() > 0.0);
    }
    #[test]
    fn test_phonon_dos_evaluate_positive() {
        let freqs: Vec<f64> = (0..50).map(|i| 1.0 + i as f64 * 0.1).collect();
        let dos = PhononDOS::from_frequencies(&freqs, 10);
        let g = dos.evaluate(2.0);
        assert!(g >= 0.0);
    }
    #[test]
    fn test_thermo_zpe_positive() {
        let freqs = vec![1e12_f64, 2e12, 3e12];
        let tp = ThermodynamicProperties::new(freqs, 300.0, 1e8);
        assert!(tp.zero_point_energy() > 0.0);
    }
    #[test]
    fn test_thermo_heat_capacity_positive_at_finite_temperature() {
        let freqs = vec![1e12_f64, 2e12, 3e12];
        let tp = ThermodynamicProperties::new(freqs, 300.0, 1e8);
        assert!(tp.heat_capacity() > 0.0);
    }
    #[test]
    fn test_thermo_heat_capacity_zero_at_zero_temperature() {
        let freqs = vec![1e12_f64, 2e12];
        let tp = ThermodynamicProperties::new(freqs, 0.0, 1e8);
        assert_eq!(tp.heat_capacity(), 0.0);
    }
    #[test]
    fn test_thermo_entropy_positive() {
        let freqs = vec![1e12_f64, 2e12];
        let tp = ThermodynamicProperties::new(freqs, 1000.0, 1e8);
        assert!(tp.entropy() > 0.0);
    }
    #[test]
    fn test_thermo_helmholtz_finite() {
        let freqs = vec![1e12_f64, 2e12];
        let tp = ThermodynamicProperties::new(freqs, 300.0, 1e8);
        assert!(tp.helmholtz_free_energy().is_finite());
    }
    #[test]
    fn test_thermo_debye_temperature_positive() {
        let freqs = vec![1e13_f64, 2e13];
        let tp = ThermodynamicProperties::new(freqs, 300.0, 0.0);
        assert!(tp.debye_temperature() > 0.0);
    }
    #[test]
    fn test_thermo_heat_capacity_dulong_petit_high_temp() {
        let n_modes = 3;
        let freqs = vec![1e10_f64; n_modes];
        let tp = ThermodynamicProperties::new(freqs, 1e6, 0.0);
        let cv = tp.heat_capacity();
        let dulong_petit = n_modes as f64 * KB;
        assert!(
            (cv - dulong_petit).abs() / dulong_petit < 0.01,
            "cv={cv:.3e}"
        );
    }
    #[test]
    fn test_gruneisen_mode_all_zero_if_no_freq_change() {
        let freqs = vec![1e12_f64, 2e12];
        let gp = GruneisenParameter::new(
            freqs.clone(),
            freqs.clone(),
            freqs.clone(),
            1e-29,
            1e-31,
            300.0,
        );
        let mg = gp.mode_gruneisen();
        assert!(mg.iter().all(|&g| g.abs() < 1e-6));
    }
    #[test]
    fn test_gruneisen_mode_positive_for_soft_mode() {
        let freqs_v = vec![1e12_f64];
        let freqs_vp = vec![0.9e12_f64];
        let freqs_vm = vec![1.1e12_f64];
        let gp = GruneisenParameter::new(freqs_v, freqs_vp, freqs_vm, 1e-29, 1e-29, 300.0);
        let mg = gp.mode_gruneisen();
        assert!(mg[0] > 0.0);
    }
    #[test]
    fn test_macroscopic_gruneisen_finite() {
        let freqs = vec![1e12_f64, 2e12];
        let freqs_p = vec![0.9e12_f64, 1.8e12];
        let freqs_m = vec![1.1e12_f64, 2.2e12];
        let gp = GruneisenParameter::new(freqs, freqs_p, freqs_m, 1e-29, 1e-31, 300.0);
        assert!(gp.macroscopic_gruneisen().is_finite());
    }
    #[test]
    fn test_gruneisen_thermal_pressure_positive() {
        let freqs = vec![1e12_f64, 2e12];
        let freqs_p = vec![0.9e12_f64, 1.8e12];
        let freqs_m = vec![1.1e12_f64, 2.2e12];
        let gp = GruneisenParameter::new(freqs, freqs_p, freqs_m, 1e-29, 1e-31, 300.0);
        let dp_dt = gp.thermal_pressure_coefficient(1e-20);
        assert!(dp_dt.is_finite());
    }
    #[test]
    fn test_qha_total_free_energies_length() {
        let ac = AnharmonicCorrections::new(
            vec![1e-29, 1.1e-29, 1.2e-29],
            vec![-1e-18, -1.01e-18, -1.005e-18],
            vec![-1e-20, -1.1e-20, -0.9e-20],
            300.0,
            1.5,
            50e9,
        );
        assert_eq!(ac.total_free_energies().len(), 3);
    }
    #[test]
    fn test_qha_equilibrium_volume_in_range() {
        let volumes: Vec<f64> = (0..5).map(|i| (1.0 + i as f64 * 0.1) * 1e-29).collect();
        let energies: Vec<f64> = volumes
            .iter()
            .map(|&v| {
                let v0 = 1.2e-29_f64;
                (v - v0).powi(2) * 1e40
            })
            .collect();
        let phonon_fe = vec![0.0_f64; 5];
        let ac = AnharmonicCorrections::new(volumes.clone(), energies, phonon_fe, 300.0, 1.5, 50e9);
        let v_eq = ac.equilibrium_volume();
        assert!(v_eq >= volumes[0] && v_eq <= volumes[volumes.len() - 1]);
    }
    #[test]
    fn test_qha_thermal_expansion_positive() {
        let ac = AnharmonicCorrections::new(
            vec![1e-29, 1.1e-29, 1.2e-29],
            vec![-1e-18, -1.1e-18, -1.05e-18],
            vec![-1e-20, -1.1e-20, -0.9e-20],
            300.0,
            1.5,
            50e9,
        );
        let cv = 3.0 * KB;
        let alpha = ac.thermal_expansion_coefficient(cv);
        assert!(alpha.is_finite());
    }
    #[test]
    fn test_phonon_frequencies_positive() {
        let dm = dynamical_matrix(4.0, 1.0, 1);
        let freqs = phonon_frequencies(&dm);
        for &w in &freqs {
            assert!(w >= 0.0);
        }
    }
    #[test]
    fn test_phonon_frequencies_value() {
        let dm = dynamical_matrix(1.0, 1.0, 1);
        let freqs = phonon_frequencies(&dm);
        assert!((freqs[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_phonon_frequencies_zero_dm() {
        let dm = [0.0_f64; 9];
        let freqs = phonon_frequencies(&dm);
        assert_eq!(freqs, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_dos_normalised() {
        let freqs = vec![1.0, 2.0, 3.0, 4.0];
        let dos = density_of_states_phonon(&freqs, 4, 5.0);
        let sum: f64 = dos.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_dos_empty() {
        let dos = density_of_states_phonon(&[], 5, 10.0);
        assert!(dos.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_dos_out_of_range_not_counted() {
        let freqs = vec![100.0];
        let dos = density_of_states_phonon(&freqs, 5, 10.0);
        assert!(dos.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_dos_bin_count() {
        let dos = density_of_states_phonon(&[1.0], 10, 10.0);
        assert_eq!(dos.len(), 10);
    }
    #[test]
    fn test_debye_temperature_formula() {
        let hbar = 1.0545718e-34_f64;
        let kb = 1.380649e-23_f64;
        let omega = 1e13_f64;
        let td = debye_temperature(omega, hbar, kb);
        let expected = hbar * omega / kb;
        assert!((td - expected).abs() < 1e-10 * expected.abs());
    }
    #[test]
    fn test_debye_temperature_positive() {
        let td = debye_temperature(1e13, 1.055e-34, 1.38e-23);
        assert!(td > 0.0);
    }
    #[test]
    fn test_debye_temperature_scales_linearly() {
        let td1 = debye_temperature(1e13, 1.0, 1.0);
        let td2 = debye_temperature(2e13, 1.0, 1.0);
        assert!((td2 / td1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_einstein_cv_high_temperature_limit() {
        let r = 8.314_f64;
        let cv = heat_capacity_einstein(10.0, 10000.0, r);
        assert!((cv - 3.0 * r).abs() < 0.01 * 3.0 * r, "cv = {cv}");
    }
    #[test]
    fn test_einstein_cv_low_temperature_zero() {
        let cv = heat_capacity_einstein(1000.0, 1.0, 8.314);
        assert!(cv < 1e-10, "cv = {cv}");
    }
    #[test]
    fn test_einstein_cv_zero_temperature() {
        let cv = heat_capacity_einstein(300.0, 0.0, 8.314);
        assert_eq!(cv, 0.0);
    }
    #[test]
    fn test_einstein_cv_positive() {
        let cv = heat_capacity_einstein(300.0, 500.0, 8.314);
        assert!(cv > 0.0);
    }
    #[test]
    fn test_einstein_cv_bounded_by_dulong_petit() {
        let r = 8.314_f64;
        let cv = heat_capacity_einstein(300.0, 500.0, r);
        assert!(cv <= 3.0 * r + 1e-10);
    }
    #[test]
    fn test_lattice_vectors_clone() {
        let lv = LatticeVectors::cubic(1.0);
        let lv2 = lv.clone();
        assert_eq!(lv, lv2);
    }
    #[test]
    fn test_unit_cell_1d_dispersion_zero_at_q0() {
        let uc = UnitCell1D::new(1.0, 1.0, 1.0);
        let omega = uc.acoustic_dispersion(0.0);
        assert!(omega.abs() < 1e-12);
    }
    #[test]
    fn test_unit_cell_1d_zone_boundary_frequency() {
        let uc = UnitCell1D::new(1.0, 4.0, 1.0);
        let wbz = uc.zone_boundary_frequency();
        assert!((wbz - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_unit_cell_1d_dispersion_at_zone_boundary() {
        let uc = UnitCell1D::new(1.0, 1.0, 1.0);
        let q_bz = PI / uc.equilibrium_spacing;
        let omega = uc.acoustic_dispersion(q_bz);
        let expected = uc.zone_boundary_frequency();
        assert!((omega - expected).abs() < 1e-10);
    }
    #[test]
    fn test_unit_cell_1d_group_velocity_at_q0() {
        let uc = UnitCell1D::new(1.0, 1.0, 1.0);
        let v = uc.group_velocity(0.0);
        let expected = uc.equilibrium_spacing * (uc.spring_const / uc.mass).sqrt();
        assert!((v - expected).abs() < 1e-12);
    }
    #[test]
    fn test_unit_cell_1d_group_velocity_zero_at_bz() {
        let uc = UnitCell1D::new(1.0, 1.0, 1.0);
        let q_bz = PI / uc.equilibrium_spacing;
        let v = uc.group_velocity(q_bz);
        assert!(v.abs() < 1e-10);
    }
    #[test]
    fn test_unit_cell_1d_dispersion_monotone() {
        let uc = UnitCell1D::new(1.0, 1.0, 1.0);
        let q_bz = PI / uc.equilibrium_spacing;
        let w1 = uc.acoustic_dispersion(0.1 * q_bz);
        let w2 = uc.acoustic_dispersion(0.5 * q_bz);
        assert!(w2 > w1);
    }
    #[test]
    fn test_diatomic_acoustic_zero_at_gamma() {
        let dc = DiatomicChain::new(1.0, 2.0, 1.0);
        let omega = dc.acoustic_frequency(0.0);
        assert!(omega.abs() < 1e-10);
    }
    #[test]
    fn test_diatomic_optical_positive() {
        let dc = DiatomicChain::new(1.0, 2.0, 1.0);
        let omega = dc.optical_frequency(0.0);
        assert!(omega > 0.0);
    }
    #[test]
    fn test_diatomic_bandgap_positive() {
        let dc = DiatomicChain::new(1.0, 2.0, 1.0);
        assert!(dc.bandgap() > 0.0);
    }
    #[test]
    fn test_diatomic_bandgap_zero_equal_masses() {
        let dc = DiatomicChain::new(1.0, 1.0, 1.0);
        assert!(dc.bandgap() < 1e-10);
    }
    #[test]
    fn test_diatomic_optical_ge_acoustic() {
        let dc = DiatomicChain::new(1.0, 3.0, 2.0);
        let q = 0.5;
        assert!(dc.optical_frequency(q) >= dc.acoustic_frequency(q));
    }
    #[test]
    fn test_debye_frequency_positive() {
        let wd = debye_frequency(1.0, 1.0, 1.0);
        assert!(wd > 0.0);
    }
    #[test]
    fn test_debye_frequency_zero_mass() {
        assert_eq!(debye_frequency(1.0, 0.0, 1.0), 0.0);
    }
    #[test]
    fn test_einstein_frequency_formula() {
        let we = einstein_frequency(4.0, 1.0);
        assert!((we - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_einstein_frequency_zero_mass() {
        assert_eq!(einstein_frequency(1.0, 0.0), 0.0);
    }
    #[test]
    fn test_phonon_dos_1d_positive() {
        let g = phonon_density_of_states_1d(0.5, 1.0);
        assert!(g > 0.0);
    }
    #[test]
    fn test_phonon_dos_1d_zero_outside() {
        assert_eq!(phonon_density_of_states_1d(1.5, 1.0), 0.0);
        assert_eq!(phonon_density_of_states_1d(-0.1, 1.0), 0.0);
    }
    #[test]
    fn test_phonon_dos_1d_diverges_near_max() {
        let g1 = phonon_density_of_states_1d(0.5, 1.0);
        let g2 = phonon_density_of_states_1d(0.99, 1.0);
        assert!(g2 > g1);
    }
    #[test]
    fn test_bose_einstein_zero_temperature() {
        let n = bose_einstein_occupation(1.0e12, 0.0, 1.055e-34);
        assert_eq!(n, 0.0);
    }
    #[test]
    fn test_bose_einstein_high_temperature() {
        let n = bose_einstein_occupation(1e10, 1e6, 1.055e-34);
        assert!(n > 0.0);
    }
    #[test]
    fn test_bose_einstein_zero_omega_infinity() {
        let n = bose_einstein_occupation(0.0, 300.0, 1.055e-34);
        assert!(n.is_infinite());
    }
    #[test]
    fn test_phonon_heat_capacity_positive() {
        let omegas = vec![1e12, 2e12, 3e12];
        let cv = phonon_heat_capacity(&omegas, 300.0, 1.055e-34);
        assert!(cv > 0.0);
    }
    #[test]
    fn test_phonon_heat_capacity_zero_temp() {
        let omegas = vec![1e12];
        let cv = phonon_heat_capacity(&omegas, 0.0, 1.055e-34);
        assert_eq!(cv, 0.0);
    }
    #[test]
    fn test_forced_oscillator_natural_frequency() {
        let fo = ForcedOscillator::new(1.0, 4.0, 0.1);
        assert!((fo.natural_frequency() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_forced_oscillator_quality_factor() {
        let fo = ForcedOscillator::new(1.0, 4.0, 0.5);
        assert!((fo.quality_factor() - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_forced_oscillator_quality_factor_zero_damping() {
        let fo = ForcedOscillator::new(1.0, 4.0, 0.0);
        assert!(fo.quality_factor().is_infinite());
    }
    #[test]
    fn test_forced_oscillator_damped_frequency_le_natural() {
        let fo = ForcedOscillator::new(1.0, 4.0, 0.5);
        assert!(fo.damped_frequency() <= fo.natural_frequency() + 1e-12);
    }
    #[test]
    fn test_forced_oscillator_steady_state_amplitude_positive() {
        let fo = ForcedOscillator::new(1.0, 4.0, 0.1);
        let amp = fo.steady_state_amplitude(1.0, 1.0);
        assert!(amp > 0.0);
    }
    #[test]
    fn test_forced_oscillator_resonance_peak() {
        let fo = ForcedOscillator::new(1.0, 4.0, 0.01);
        let w0 = fo.natural_frequency();
        let amp_res = fo.steady_state_amplitude(1.0, w0);
        let amp_off = fo.steady_state_amplitude(1.0, 0.1 * w0);
        assert!(amp_res > amp_off);
    }
    #[test]
    fn test_gruneisen_parameter_formula() {
        let g = gruneisen_parameter(2.0, 3.0, 6.0, 1.0);
        assert!((g - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_gruneisen_parameter_zero_cv() {
        let g = gruneisen_parameter(1.0, 1.0, 0.0, 1.0);
        assert_eq!(g, 0.0);
    }
}
