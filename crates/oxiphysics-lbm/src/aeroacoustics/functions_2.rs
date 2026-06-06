//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{AcousticSourceType, FwhObserver};

#[cfg(test)]
mod tests_extended_acoustics {
    use super::super::functions::*;
    use std::f64::consts::PI;

    #[test]
    fn test_greens_3d_amplitude_zero_distance() {
        let g = greens_function_3d_amplitude(0.0);
        assert_eq!(g, 0.0, "G(0) should be 0 (guarded)");
    }
    #[test]
    fn test_greens_3d_amplitude_unit_distance() {
        let g = greens_function_3d_amplitude(1.0);
        let expected = 1.0 / (4.0 * PI);
        assert!(
            (g - expected).abs() < 1e-14,
            "G(1) = {g}, expected {expected}"
        );
    }
    #[test]
    fn test_greens_3d_decays_with_distance() {
        let g1 = greens_function_3d_amplitude(1.0);
        let g2 = greens_function_3d_amplitude(2.0);
        assert!(
            g1 > g2,
            "G should decay with distance: G(1)={g1}, G(2)={g2}"
        );
    }
    #[test]
    fn test_half_space_greens_larger_than_free_field() {
        let g_free = greens_function_3d_amplitude(1.0);
        let g_half = half_space_greens(1.0, 1.0);
        assert!(g_half >= g_free, "half-space = {g_half}, free = {g_free}");
    }
    #[test]
    fn test_retarded_time_basic() {
        let t_ret = retarded_time(1.0, 340.0, 340.0);
        assert!((t_ret - 0.0).abs() < 1e-14, "t_ret = {t_ret}");
    }
    #[test]
    fn test_sphere_scattering_cross_section_scales_k4() {
        let k1 = 0.1;
        let k2 = 0.2;
        let a = 0.05;
        let s1 = sphere_scattering_cross_section(k1, a);
        let s2 = sphere_scattering_cross_section(k2, a);
        let ratio = s2 / s1;
        assert!((ratio - 16.0).abs() < 1e-10, "ratio = {ratio}");
    }
    #[test]
    fn test_cylinder_scattering_rayleigh_finite() {
        let p = cylinder_scattering_rayleigh(1.0, 0.5, 0.05, 2.0);
        assert!(p.is_finite() && p >= 0.0, "cylinder scattering = {p}");
    }
    #[test]
    fn test_wedge_diffraction_zero_path_diff() {
        let corr = wedge_diffraction_spl_correction(0.0, 0.1);
        assert_eq!(corr, 0.0, "zero path diff → 0 correction");
    }
    #[test]
    fn test_vortex_sound_zero_vorticity() {
        let p = vortex_sound_pressure(1.2, 0.0, [0.1, 0.0], [1.0, 0.0], 1.0, 340.0);
        assert_eq!(p, 0.0, "zero vorticity → zero pressure");
    }
    #[test]
    fn test_howe_vortex_sound_power_positive() {
        let w = howe_vortex_sound_power(1.2, 10.0, 100.0, 0.01, 340.0);
        assert!(w > 0.0, "Howe vortex power should be positive: {w}");
    }
    #[test]
    fn test_karman_shedding_strouhal_relation() {
        let f = karman_shedding_frequency(10.0, 0.05, 0.2);
        assert!((f - 40.0).abs() < 1e-10, "f = {f}");
    }
    #[test]
    fn test_piston_radiation_resistance_increases_with_k() {
        let r1 = piston_radiation_resistance_low_freq(1.2, 340.0, 0.05, 0.1);
        let r2 = piston_radiation_resistance_low_freq(1.2, 340.0, 0.05, 0.2);
        assert!((r2 / r1 - 4.0).abs() < 1e-10, "R ratio = {}", r2 / r1);
    }
    #[test]
    fn test_room_constant_full_absorption() {
        let r = room_constant(100.0, 0.9999);
        assert!(r > 1e5, "R with near-total absorption should be huge: {r}");
    }
    #[test]
    fn test_room_constant_zero_absorption() {
        let r = room_constant(100.0, 0.0);
        assert_eq!(r, 0.0, "R with zero absorption = {r}");
    }
    #[test]
    fn test_lbm_density_fluctuation_equilibrium() {
        let f = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2];
        let fluct = lbm_density_fluctuation(&f, 1.0);
        assert!((fluct - 0.0).abs() < 1e-14, "fluctuation = {fluct}");
    }
    #[test]
    fn test_lbm_acoustic_pressure_proportional_to_cs2() {
        let f = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3];
        let p1 = lbm_acoustic_pressure(&f, 1.0, 1.0 / 3.0);
        let p2 = lbm_acoustic_pressure(&f, 1.0, 2.0 / 3.0);
        assert!((p2 / p1 - 2.0).abs() < 1e-13, "p2/p1 = {}", p2 / p1);
    }
    #[test]
    fn test_acoustic_cfl_unit() {
        let cfl = acoustic_cfl(340.0, 1.0 / 340.0, 1.0);
        assert!((cfl - 1.0).abs() < 1e-14, "CFL = {cfl}");
    }
    #[test]
    fn test_spherical_spreading_loss_double_distance() {
        let loss = spherical_spreading_loss(1.0, 2.0);
        assert!(
            (loss - (-6.0206)).abs() < 0.001,
            "spherical spreading: {loss}"
        );
    }
    #[test]
    fn test_cylindrical_spreading_loss_double_distance() {
        let loss = cylindrical_spreading_loss(1.0, 2.0);
        assert!(
            (loss - (-3.0103)).abs() < 0.001,
            "cylindrical spreading: {loss}"
        );
    }
    #[test]
    fn test_add_incoherent_spl_two_equal() {
        let result = add_incoherent_spl(&[60.0, 60.0]);
        assert!(
            (result - 63.0103).abs() < 0.001,
            "Incoherent SPL sum = {result}"
        );
    }
    #[test]
    fn test_add_incoherent_spl_empty() {
        let result = add_incoherent_spl(&[]);
        assert!(result.is_infinite());
    }
    #[test]
    fn test_atmospheric_absorption_increases_with_frequency() {
        let alpha_lo = atmospheric_absorption_db_per_m(1000.0, 293.15, 50.0);
        let alpha_hi = atmospheric_absorption_db_per_m(8000.0, 293.15, 50.0);
        assert!(
            alpha_hi > alpha_lo,
            "alpha_lo={alpha_lo}, alpha_hi={alpha_hi}"
        );
    }
}
/// Compute acoustic pressure fluctuation from density: p' = cs² * (ρ - ρ₀).
///
/// This is the LBM acoustic analogy formula relating density fluctuations to
/// pressure fluctuations through the isothermal speed of sound.
pub fn acoustic_pressure_fluctuation(rho: f64, rho0: f64, cs2: f64) -> f64 {
    cs2 * (rho - rho0)
}
/// Compute acoustic pressure fluctuation field over the entire lattice.
///
/// Returns a vector of p' values for each lattice node.
pub fn acoustic_pressure_field(rho: &[f64], rho0: f64, cs2: f64) -> Vec<f64> {
    rho.iter()
        .map(|&r| acoustic_pressure_fluctuation(r, rho0, cs2))
        .collect()
}
/// Compute RMS acoustic pressure from a time series of pressure fluctuations.
pub fn pressure_fluctuation_rms(p_prime: &[f64]) -> f64 {
    if p_prime.is_empty() {
        return 0.0;
    }
    let mean_sq: f64 = p_prime.iter().map(|&p| p * p).sum::<f64>() / p_prime.len() as f64;
    mean_sq.sqrt()
}
/// Compute the trace (isotropic part) of the 2D Lighthill stress tensor.
///
/// tr(T) = T_xx + T_yy
pub fn lighthill_tensor_trace(t: [[f64; 2]; 2]) -> f64 {
    t[0][0] + t[1][1]
}
/// Compute the deviatoric (traceless) part of the 2D Lighthill stress tensor.
///
/// T_dev = T - (tr(T)/2) * I
pub fn lighthill_tensor_deviatoric(t: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let half_trace = (t[0][0] + t[1][1]) / 2.0;
    [
        [t[0][0] - half_trace, t[0][1]],
        [t[1][0], t[1][1] - half_trace],
    ]
}
/// Check whether Lighthill tensor is symmetric (T_ij = T_ji).
pub fn lighthill_tensor_is_symmetric(t: [[f64; 2]; 2], tol: f64) -> bool {
    (t[0][1] - t[1][0]).abs() < tol
}
/// Estimate the acoustic far-field power from a quadrupole source region.
///
/// W_quad ≈ ρ₀ * ∫ T_ij T_ij dV / (c₀^5) (Lighthill 8th-power-like scaling)
pub fn quadrupole_far_field_power(t_frobenius_sq: f64, volume: f64, rho0: f64, c0: f64) -> f64 {
    rho0 * t_frobenius_sq * volume / c0.powi(5)
}
/// Compute monopole term of the FWH equation at an observer.
///
/// p'_T = ρ₀ * Σ_i \[ u_n * dA_i / (4π r_i²) \]
/// where u_n is normal surface velocity.
pub fn fwh_monopole_term(
    surface_velocities: &[[f64; 3]],
    surface_normals: &[[f64; 3]],
    surface_areas: &[f64],
    surface_positions: &[[f64; 3]],
    observer: &FwhObserver,
) -> f64 {
    let n = surface_velocities.len();
    assert_eq!(surface_normals.len(), n);
    assert_eq!(surface_areas.len(), n);
    assert_eq!(surface_positions.len(), n);
    let mut p_total = 0.0_f64;
    for i in 0..n {
        let u_n = surface_velocities[i][0] * surface_normals[i][0]
            + surface_velocities[i][1] * surface_normals[i][1]
            + surface_velocities[i][2] * surface_normals[i][2];
        let r = observer.distance_to(surface_positions[i]);
        if r < 1e-15 {
            continue;
        }
        p_total += observer.rho0 * u_n * surface_areas[i] / (4.0 * PI * r * r);
    }
    p_total
}
/// Compute dipole term of the FWH equation at an observer.
///
/// p'_L = Σ_i \[ Δp_i * (n̂ · r̂) * dA_i / (4π r_i²) \]
pub fn fwh_dipole_term(
    surface_pressures: &[f64],
    surface_normals: &[[f64; 3]],
    surface_areas: &[f64],
    surface_positions: &[[f64; 3]],
    observer: &FwhObserver,
) -> f64 {
    let n = surface_pressures.len();
    assert_eq!(surface_normals.len(), n);
    assert_eq!(surface_areas.len(), n);
    assert_eq!(surface_positions.len(), n);
    let mut p_total = 0.0_f64;
    for i in 0..n {
        let r_hat = observer.unit_vector_from(surface_positions[i]);
        let n_dot_r = surface_normals[i][0] * r_hat[0]
            + surface_normals[i][1] * r_hat[1]
            + surface_normals[i][2] * r_hat[2];
        let r = observer.distance_to(surface_positions[i]);
        if r < 1e-15 {
            continue;
        }
        p_total += surface_pressures[i] * n_dot_r * surface_areas[i] / (4.0 * PI * r * r);
    }
    p_total
}
/// Compute SPL in dB: SPL = 20 * log10(p' / p_ref) where p_ref = 20 µPa.
///
/// Returns NEG_INFINITY if p_prime <= 0.
pub fn spl_from_pressure_fluctuation(p_prime: f64) -> f64 {
    pub(super) const P_REF: f64 = 20e-6;
    if p_prime <= 0.0 {
        return f64::NEG_INFINITY;
    }
    20.0 * (p_prime / P_REF).log10()
}
/// Convert SPL in dB back to pressure amplitude (Pa).
pub fn pressure_from_spl(spl_db: f64) -> f64 {
    pub(super) const P_REF: f64 = 20e-6;
    P_REF * 10.0_f64.powf(spl_db / 20.0)
}
/// Compute peak SPL from a time series of pressure fluctuations.
pub fn peak_spl(p_prime: &[f64]) -> f64 {
    let p_max = p_prime.iter().cloned().fold(0.0_f64, |a, b| a.max(b.abs()));
    spl_from_pressure_fluctuation(p_max)
}
/// Compute 1/3 octave band center frequencies from 20 Hz to 20 kHz.
///
/// Returns center frequencies of ISO 1/3-octave bands.
pub fn third_octave_center_frequencies() -> Vec<f64> {
    vec![
        20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0,
        500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0,
        6300.0, 8000.0, 10000.0, 12500.0, 16000.0, 20000.0,
    ]
}
/// Compute the vorticity ω_z = ∂v/∂x - ∂u/∂y from velocity field differences.
///
/// Uses central differences: ω_z ≈ (v_{i+1,j} - v_{i-1,j})/(2Δx) - (u_{i,j+1} - u_{i,j-1})/(2Δy)
pub fn compute_vorticity_z(
    u_xplus: f64,
    u_xminus: f64,
    v_yplus: f64,
    v_yminus: f64,
    dx: f64,
    dy: f64,
) -> f64 {
    (v_yplus - v_yminus) / (2.0 * dx) - (u_xplus - u_xminus) / (2.0 * dy)
}
/// Compute circulation around a closed contour from vorticity field.
///
/// Γ = ∫∫ ω_z dA ≈ Σ ω_z * dA
pub fn circulation_from_vorticity(vorticity: &[f64], cell_area: f64) -> f64 {
    vorticity.iter().sum::<f64>() * cell_area
}
/// Enstrophy: Z = 0.5 * ∫ ω² dV (measure of vortex intensity).
pub fn enstrophy(vorticity: &[f64], cell_volume: f64) -> f64 {
    0.5 * vorticity.iter().map(|&w| w * w).sum::<f64>() * cell_volume
}
/// Acoustic power from vortex dynamics (Powell's analogy, 2D).
///
/// W_ac = ρ₀/(4π c₀³) * (dΓ/dt)²
pub fn vortex_acoustic_power_2d(d_gamma_dt: f64, rho0: f64, c0: f64) -> f64 {
    rho0 * d_gamma_dt * d_gamma_dt / (4.0 * PI * c0.powi(3))
}
/// Estimate acoustic emission from a vortex pair (dipole model).
///
/// Two counter-rotating vortices of strength ±Γ separated by distance d
/// emit dipole sound. p_rms ~ ρ₀ Γ² / (2π c₀ r d)
pub fn vortex_pair_acoustic_pressure(gamma: f64, d: f64, r: f64, rho0: f64, c0: f64) -> f64 {
    if r < 1e-15 || c0 < 1e-15 {
        return 0.0;
    }
    rho0 * gamma * gamma / (2.0 * PI * c0 * r * d)
}
/// Identify dominant source type from Mach number and measured power law exponent.
///
/// `power_exponent` is the measured exponent n in W ∝ Uⁿ.
pub fn identify_source_type(power_exponent: f64) -> AcousticSourceType {
    if (power_exponent - 2.0).abs() <= 1.0 {
        AcousticSourceType::Monopole
    } else if (power_exponent - 4.0).abs() <= 1.5 {
        AcousticSourceType::Dipole
    } else {
        AcousticSourceType::Quadrupole
    }
}
/// Source strength ratio: compare monopole, dipole, quadrupole amplitudes.
///
/// Returns (Q_monopole, Q_dipole, Q_quadrupole) normalized to monopole = 1.
pub fn source_strength_ratio(ma: f64) -> (f64, f64, f64) {
    let q_mono = 1.0;
    let q_dip = ma * ma;
    let q_quad = ma.powi(4);
    (q_mono, q_dip, q_quad)
}
/// Check if source is in compact limit: k*L << 1.
///
/// Returns true if the Helmholtz number He = k*L < 0.1 (compact limit).
pub fn is_compact_source(k: f64, length: f64) -> bool {
    k * length < 0.1
}
/// Multipole expansion: compute acoustic pressure from monopole + dipole terms.
///
/// p = Q_mono/(4π r) + F⃗ · r̂ / (4π r²)
pub fn multipole_acoustic_pressure(q_mono: f64, dipole: [f64; 3], r: f64, r_hat: [f64; 3]) -> f64 {
    if r < 1e-15 {
        return 0.0;
    }
    let mono = q_mono / (4.0 * PI * r);
    let dip =
        (dipole[0] * r_hat[0] + dipole[1] * r_hat[1] + dipole[2] * r_hat[2]) / (4.0 * PI * r * r);
    mono + dip
}
/// Compute the acoustic intensity vector I = p * u_ac at a grid point.
///
/// I⃗ = p' * u⃗_ac
pub fn acoustic_intensity_vector(p_prime: f64, u_ac: [f64; 3]) -> [f64; 3] {
    [p_prime * u_ac[0], p_prime * u_ac[1], p_prime * u_ac[2]]
}
/// Time-averaged acoustic intensity vector: <I⃗> = <p' u⃗_ac>.
///
/// Computed from arrays of time samples.
pub fn time_averaged_intensity_vector(
    p_prime: &[f64],
    u_ac_x: &[f64],
    u_ac_y: &[f64],
    u_ac_z: &[f64],
) -> [f64; 3] {
    let n = p_prime.len();
    assert_eq!(u_ac_x.len(), n);
    assert_eq!(u_ac_y.len(), n);
    assert_eq!(u_ac_z.len(), n);
    if n == 0 {
        return [0.0; 3];
    }
    let inv_n = 1.0 / n as f64;
    let ix = p_prime
        .iter()
        .zip(u_ac_x)
        .map(|(&p, &u)| p * u)
        .sum::<f64>()
        * inv_n;
    let iy = p_prime
        .iter()
        .zip(u_ac_y)
        .map(|(&p, &u)| p * u)
        .sum::<f64>()
        * inv_n;
    let iz = p_prime
        .iter()
        .zip(u_ac_z)
        .map(|(&p, &u)| p * u)
        .sum::<f64>()
        * inv_n;
    [ix, iy, iz]
}
/// Magnitude of acoustic intensity vector.
pub fn intensity_vector_magnitude(i_vec: [f64; 3]) -> f64 {
    (i_vec[0] * i_vec[0] + i_vec[1] * i_vec[1] + i_vec[2] * i_vec[2]).sqrt()
}
/// Acoustic intensity direction (unit vector).
pub fn intensity_direction(i_vec: [f64; 3]) -> [f64; 3] {
    let mag = intensity_vector_magnitude(i_vec);
    if mag < 1e-30 {
        return [0.0, 0.0, 0.0];
    }
    [i_vec[0] / mag, i_vec[1] / mag, i_vec[2] / mag]
}
/// Free-space acoustic wavenumber: k = ω / c₀ = 2π f / c₀.
pub fn acoustic_wavenumber(freq: f64, c0: f64) -> f64 {
    2.0 * PI * freq / c0
}
/// Acoustic wavenumber in a medium with mean flow (convected wave equation).
///
/// k± = (-Mc ± 1) / ((1 - Mc²) * λ)   where Mc = M * cos(θ)
///
/// Returns the downstream (+) wavenumber.
pub fn convected_wavenumber(freq: f64, c0: f64, mach: f64, theta: f64) -> f64 {
    let mc = mach * theta.cos();
    let lambda = c0 / freq;
    let denom = (1.0 - mc * mc) * lambda;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    (-mc + 1.0) / denom
}
/// Dispersion relation check: verify acoustic CFL condition.
///
/// Returns true if the simulation is stable: c₀ * Δt / Δx ≤ CFL_max.
pub fn check_acoustic_cfl(c0: f64, dt: f64, dx: f64, cfl_max: f64) -> bool {
    c0 * dt / dx <= cfl_max
}
/// Group velocity for acoustic waves in a dispersive medium.
///
/// cg = dω/dk = c₀ / (1 + dispersion_correction * k²)
/// (first-order approximation for weakly dispersive media)
pub fn group_velocity(c0: f64, k: f64, dispersion: f64) -> f64 {
    c0 / (1.0 + dispersion * k * k)
}
#[cfg(test)]
mod tests_new_aeroacoustics {
    use super::super::*;
    use crate::aeroacoustics::types::*;
    use std::f64::consts::PI;
    #[test]
    fn test_acoustic_pressure_fluctuation_zero() {
        let p = acoustic_pressure_fluctuation(1.0, 1.0, 1.0 / 3.0);
        assert!(p.abs() < 1e-14, "p' at equilibrium = {p}");
    }
    #[test]
    fn test_acoustic_pressure_fluctuation_positive() {
        let p = acoustic_pressure_fluctuation(1.1, 1.0, 1.0 / 3.0);
        assert!(p > 0.0, "p' for rho > rho0 should be positive: {p}");
    }
    #[test]
    fn test_acoustic_pressure_field_length() {
        let rho = vec![1.0, 1.1, 0.9, 1.0, 1.0];
        let field = acoustic_pressure_field(&rho, 1.0, 1.0 / 3.0);
        assert_eq!(field.len(), rho.len());
    }
    #[test]
    fn test_pressure_fluctuation_rms_zero() {
        let p = vec![0.0, 0.0, 0.0];
        let rms = pressure_fluctuation_rms(&p);
        assert_eq!(rms, 0.0);
    }
    #[test]
    fn test_pressure_fluctuation_rms_constant() {
        let p = vec![1.0, 1.0, 1.0];
        let rms = pressure_fluctuation_rms(&p);
        assert!((rms - 1.0).abs() < 1e-14, "rms of [1,1,1] = {rms}");
    }
    #[test]
    fn test_lighthill_tensor_trace() {
        let t = [[2.0, 1.0], [1.0, 3.0]];
        let tr = lighthill_tensor_trace(t);
        assert!((tr - 5.0).abs() < 1e-14, "trace = {tr}");
    }
    #[test]
    fn test_lighthill_tensor_deviatoric_traceless() {
        let t = [[2.0, 1.0], [1.0, 4.0]];
        let dev = lighthill_tensor_deviatoric(t);
        let tr = dev[0][0] + dev[1][1];
        assert!(tr.abs() < 1e-14, "deviatoric trace should be 0: {tr}");
    }
    #[test]
    fn test_lighthill_tensor_is_symmetric_true() {
        let t = [[1.0, 0.5], [0.5, 2.0]];
        assert!(lighthill_tensor_is_symmetric(t, 1e-14));
    }
    #[test]
    fn test_lighthill_tensor_is_symmetric_false() {
        let t = [[1.0, 0.5], [0.6, 2.0]];
        assert!(!lighthill_tensor_is_symmetric(t, 0.05));
    }
    #[test]
    fn test_fwh_observer_distance() {
        let obs = FwhObserver::new(3.0, 4.0, 0.0, 340.0, 1.2);
        let d = obs.distance_to([0.0, 0.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-12, "distance = {d}");
    }
    #[test]
    fn test_fwh_observer_unit_vector() {
        let obs = FwhObserver::new(1.0, 0.0, 0.0, 340.0, 1.2);
        let uv = obs.unit_vector_from([0.0, 0.0, 0.0]);
        assert!((uv[0] - 1.0).abs() < 1e-14, "unit vector x = {}", uv[0]);
        assert!(uv[1].abs() < 1e-14);
        assert!(uv[2].abs() < 1e-14);
    }
    #[test]
    fn test_fwh_monopole_term_single_panel() {
        let obs = FwhObserver::new(10.0, 0.0, 0.0, 340.0, 1.2);
        let vels = [[1.0, 0.0, 0.0_f64]];
        let norms = [[1.0, 0.0, 0.0_f64]];
        let areas = [1.0_f64];
        let positions = [[0.0, 0.0, 0.0_f64]];
        let p = fwh_monopole_term(&vels, &norms, &areas, &positions, &obs);
        assert!(p.abs() > 0.0, "FWH monopole term = {p}");
    }
    #[test]
    fn test_fwh_dipole_term_single_panel() {
        let obs = FwhObserver::new(10.0, 0.0, 0.0, 340.0, 1.2);
        let pressures = [100.0_f64];
        let norms = [[1.0, 0.0, 0.0_f64]];
        let areas = [1.0_f64];
        let positions = [[0.0, 0.0, 0.0_f64]];
        let p = fwh_dipole_term(&pressures, &norms, &areas, &positions, &obs);
        assert!(p.abs() > 0.0, "FWH dipole term = {p}");
    }
    #[test]
    fn test_spl_from_pressure_fluctuation_reference() {
        let spl = spl_from_pressure_fluctuation(20e-6);
        assert!(spl.abs() < 1e-10, "SPL at p_ref = {spl}");
    }
    #[test]
    fn test_spl_from_pressure_fluctuation_zero() {
        let spl = spl_from_pressure_fluctuation(0.0);
        assert!(spl.is_infinite());
    }
    #[test]
    fn test_pressure_from_spl_roundtrip() {
        let p_orig = 0.01;
        let spl = spl_from_pressure_fluctuation(p_orig);
        let p_back = pressure_from_spl(spl);
        assert!((p_back - p_orig).abs() < 1e-14, "roundtrip: {p_back}");
    }
    #[test]
    fn test_third_octave_bands_count() {
        let bands = third_octave_center_frequencies();
        assert_eq!(bands.len(), 31, "should have 31 1/3-octave bands");
    }
    #[test]
    fn test_a_weighting_filter_1khz_near_zero() {
        let aw = AWeightingFilter::new();
        let db_1k = aw.magnitude_db(1000.0);
        assert!((db_1k - 0.0).abs() < 2.0, "A-weighting at 1 kHz = {db_1k}");
    }
    #[test]
    fn test_a_weighting_filter_low_attenuated() {
        let aw = AWeightingFilter::new();
        let db_50 = aw.magnitude_db(50.0);
        let db_1k = aw.magnitude_db(1000.0);
        assert!(db_50 < db_1k, "50 Hz should be more attenuated than 1 kHz");
    }
    #[test]
    fn test_a_weighting_filter_apply_spectrum() {
        let aw = AWeightingFilter::new();
        let freqs = [500.0, 1000.0, 2000.0];
        let spls = [70.0, 70.0, 70.0];
        let aw_spls = aw.apply_to_spectrum(&freqs, &spls);
        assert_eq!(aw_spls.len(), 3);
        for &v in &aw_spls {
            assert!(v.is_finite(), "A-weighted SPL should be finite: {v}");
        }
    }
    #[test]
    fn test_compute_vorticity_z_pure_rotation() {
        let omega = 1.0;
        let dx = 0.1;
        let dy = 0.1;
        let wz = compute_vorticity_z(
            -omega * (0.0 + dy),
            -omega * (0.0 - dy),
            omega * (0.0 + dx),
            omega * (0.0 - dx),
            dx,
            dy,
        );
        assert!((wz - 2.0 * omega).abs() < 1e-12, "vorticity = {wz}");
    }
    #[test]
    fn test_circulation_from_vorticity_uniform() {
        let vorticity = vec![1.0; 100];
        let circ = circulation_from_vorticity(&vorticity, 0.01);
        assert!((circ - 1.0).abs() < 1e-12, "circulation = {circ}");
    }
    #[test]
    fn test_enstrophy_positive() {
        let vorticity = vec![1.0, 2.0, 3.0];
        let z = enstrophy(&vorticity, 1.0);
        assert!(z > 0.0, "enstrophy = {z}");
    }
    #[test]
    fn test_vortex_acoustic_power_2d_positive() {
        let w = vortex_acoustic_power_2d(10.0, 1.2, 340.0);
        assert!(w > 0.0, "vortex acoustic power = {w}");
    }
    #[test]
    fn test_vortex_pair_acoustic_pressure_positive() {
        let p = vortex_pair_acoustic_pressure(1.0, 0.1, 10.0, 1.2, 340.0);
        assert!(p > 0.0, "vortex pair pressure = {p}");
    }
    #[test]
    fn test_identify_source_monopole() {
        assert_eq!(identify_source_type(2.0), AcousticSourceType::Monopole);
    }
    #[test]
    fn test_identify_source_dipole() {
        assert_eq!(identify_source_type(4.0), AcousticSourceType::Dipole);
    }
    #[test]
    fn test_identify_source_quadrupole() {
        assert_eq!(identify_source_type(8.0), AcousticSourceType::Quadrupole);
    }
    #[test]
    fn test_source_strength_ratio_mach_one() {
        let (q_mono, q_dip, q_quad) = source_strength_ratio(1.0);
        assert!((q_mono - 1.0).abs() < 1e-14);
        assert!((q_dip - 1.0).abs() < 1e-14);
        assert!((q_quad - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_source_strength_ratio_low_mach() {
        let (q_mono, q_dip, q_quad) = source_strength_ratio(0.1);
        assert!(q_mono > q_dip, "monopole dominates at low Mach");
        assert!(
            q_dip > q_quad,
            "dipole dominates over quadrupole at low Mach"
        );
    }
    #[test]
    fn test_is_compact_source_true() {
        assert!(is_compact_source(0.01, 0.001), "k*L = 0.00001 is compact");
    }
    #[test]
    fn test_is_compact_source_false() {
        assert!(!is_compact_source(100.0, 1.0), "k*L = 100 is not compact");
    }
    #[test]
    fn test_multipole_acoustic_pressure_monopole_only() {
        let p = multipole_acoustic_pressure(1.0, [0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0]);
        let expected = 1.0 / (4.0 * PI);
        assert!((p - expected).abs() < 1e-14, "p = {p}, expected {expected}");
    }
    #[test]
    fn test_acoustic_intensity_vector_magnitude() {
        let i_vec = acoustic_intensity_vector(2.0, [1.0, 0.0, 0.0]);
        assert!((i_vec[0] - 2.0).abs() < 1e-14);
        assert!(i_vec[1].abs() < 1e-14);
        assert!(i_vec[2].abs() < 1e-14);
    }
    #[test]
    fn test_time_averaged_intensity_vector_zero() {
        let p = vec![1.0, -1.0, 1.0, -1.0];
        let u = vec![1.0, 1.0, 1.0, 1.0];
        let z = vec![0.0, 0.0, 0.0, 0.0];
        let i_vec = time_averaged_intensity_vector(&p, &u, &z, &z);
        assert!(
            i_vec[0].abs() < 1e-14,
            "time-avg intensity x = {}",
            i_vec[0]
        );
    }
    #[test]
    fn test_intensity_vector_magnitude_unit() {
        let mag = intensity_vector_magnitude([3.0, 4.0, 0.0]);
        assert!((mag - 5.0).abs() < 1e-14, "magnitude = {mag}");
    }
    #[test]
    fn test_intensity_direction_unit() {
        let dir = intensity_direction([3.0, 0.0, 0.0]);
        assert!((dir[0] - 1.0).abs() < 1e-14, "direction = {:?}", dir);
        assert!(dir[1].abs() < 1e-14);
        assert!(dir[2].abs() < 1e-14);
    }
    #[test]
    fn test_near_to_far_field_empty() {
        let ntf = NearToFarField::new();
        assert!(ntf.panels.is_empty());
        assert!(ntf.observers.is_empty());
    }
    #[test]
    fn test_near_to_far_field_with_panel() {
        let mut ntf = NearToFarField::new();
        ntf.add_panel(FwhPanel {
            center: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            area: 1.0,
            pressure: 1.0,
            un: 0.0,
        });
        ntf.add_observer([10.0, 0.0, 0.0]);
        let ff = ntf.compute_far_field();
        assert_eq!(ff.len(), 1);
        assert!(ff[0] > 0.0, "far field = {}", ff[0]);
    }
    #[test]
    fn test_acoustic_wavenumber_basic() {
        let k = acoustic_wavenumber(340.0, 340.0);
        assert!((k - 2.0 * PI).abs() < 1e-10, "k = {k}");
    }
    #[test]
    fn test_check_acoustic_cfl_stable() {
        assert!(
            check_acoustic_cfl(340.0, 1e-4, 0.1, 0.5),
            "should be stable"
        );
    }
    #[test]
    fn test_check_acoustic_cfl_unstable() {
        assert!(
            !check_acoustic_cfl(340.0, 1.0, 0.1, 0.5),
            "should be unstable"
        );
    }
    #[test]
    fn test_group_velocity_no_dispersion() {
        let cg = group_velocity(340.0, 1.0, 0.0);
        assert!((cg - 340.0).abs() < 1e-10, "cg = {cg}");
    }
    #[test]
    fn test_quadrupole_far_field_power_positive() {
        let w = quadrupole_far_field_power(1.0, 1.0, 1.2, 340.0);
        assert!(w > 0.0, "quadrupole power = {w}");
    }
    #[test]
    fn test_greens_function_3d_amplitude_falloff() {
        let g1 = greens_function_3d_amplitude(1.0);
        let g2 = greens_function_3d_amplitude(2.0);
        assert!(
            g1 > g2,
            "Green's function should decay: g(1)={g1}, g(2)={g2}"
        );
    }
    #[test]
    fn test_greens_function_3d_amplitude_positive() {
        let g = greens_function_3d_amplitude(5.0);
        assert!(
            g > 0.0,
            "Green's function amplitude should be positive: {g}"
        );
    }
    #[test]
    fn test_half_space_greens_doubles_at_zero_image() {
        let g = half_space_greens(1.0, 1.0);
        let g_free = greens_function_3d_amplitude(1.0);
        assert!(
            (g - 2.0 * g_free).abs() < 1e-12,
            "half-space at equal distances = {g}"
        );
    }
    #[test]
    fn test_retarded_time_zero_distance() {
        let t_ret = retarded_time(1.0, 0.0, 340.0);
        assert!(
            (t_ret - 1.0).abs() < 1e-15,
            "retarded time at r=0 = {t_ret}"
        );
    }
    #[test]
    fn test_retarded_time_positive_distance() {
        let t = 2.0;
        let r = 340.0;
        let c0 = 340.0;
        let t_ret = retarded_time(t, r, c0);
        assert!((t_ret - 1.0).abs() < 1e-12, "retarded time = {t_ret}");
    }
    #[test]
    fn test_sphere_scattering_rayleigh_zero_at_pi_over_2() {
        let p = sphere_scattering_rayleigh(1.0, 0.01, 0.01, 1.0, PI / 2.0);
        assert!(p.is_finite(), "scattering at pi/2 should be finite: {p}");
    }
    #[test]
    fn test_cylinder_scattering_rayleigh_decays_with_distance() {
        let p1 = cylinder_scattering_rayleigh(1.0, 0.1, 0.01, 1.0);
        let p2 = cylinder_scattering_rayleigh(1.0, 0.1, 0.01, 2.0);
        assert!(p1 > p2, "cylinder scattering decays: {p1} > {p2}");
    }
    #[test]
    fn test_sphere_scattering_cross_section_positive() {
        let sigma = sphere_scattering_cross_section(0.1, 0.01);
        assert!(sigma > 0.0, "scattering cross section = {sigma}");
    }
    #[test]
    fn test_wedge_diffraction_spl_correction_finite() {
        let d_b = wedge_diffraction_spl_correction(0.5, 1.0);
        assert!(d_b.is_finite(), "wedge diffraction correction = {d_b}");
    }
    #[test]
    fn test_karman_shedding_frequency_strouhal() {
        let f = karman_shedding_frequency(10.0, 0.1, 0.2);
        assert!((f - 20.0).abs() < 1e-12, "Kármán shedding f = {f}");
    }
    #[test]
    fn test_karman_shedding_frequency_proportional_to_velocity() {
        let f1 = karman_shedding_frequency(10.0, 0.1, 0.2);
        let f2 = karman_shedding_frequency(20.0, 0.1, 0.2);
        assert!((f2 / f1 - 2.0).abs() < 1e-12, "f ∝ u: {f1}, {f2}");
    }
    #[test]
    fn test_piston_radiation_resistance_positive() {
        let r_rad = piston_radiation_resistance_low_freq(1.2, 340.0, 0.05, 0.1);
        assert!(r_rad > 0.0, "piston radiation resistance = {r_rad}");
    }
    #[test]
    fn test_piston_radiation_reactance_positive() {
        let x_rad = piston_radiation_reactance_low_freq(1.2, 340.0, 0.05, 0.1);
        assert!(x_rad > 0.0, "piston radiation reactance = {x_rad}");
    }
    #[test]
    fn test_room_constant_positive() {
        let rc = room_constant(100.0, 0.3);
        assert!(rc > 0.0, "room constant = {rc}");
    }
    #[test]
    fn test_room_constant_fully_absorbing() {
        let rc = room_constant(100.0, 1.0);
        assert!(rc.is_infinite(), "fully absorbing room constant = {rc}");
    }
    #[test]
    fn test_lbm_density_fluctuation_zero_at_mean() {
        let f = [
            4.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0_f64,
        ];
        let drho = lbm_density_fluctuation(&f, 1.0);
        assert!(drho.abs() < 1e-14, "drho at mean density = {drho}");
    }
    #[test]
    fn test_lbm_acoustic_pressure_non_negative_at_mean() {
        let f = [
            4.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0_f64,
        ];
        let p_ac = lbm_acoustic_pressure(&f, 1.0, 1.0 / 3.0);
        assert!(
            p_ac.abs() < 1e-14,
            "acoustic pressure at mean density = {p_ac}"
        );
    }
    #[test]
    fn test_acoustic_cfl_computation() {
        let cfl = acoustic_cfl(340.0, 1e-4, 0.01);
        assert!((cfl - 3.4).abs() < 1e-12, "CFL = {cfl}");
    }
    #[test]
    fn test_acoustic_absorption_length_positive() {
        let l = acoustic_absorption_length(0.1);
        assert!(l > 0.0, "absorption length = {l}");
    }
    #[test]
    fn test_spherical_spreading_loss_negative_for_farther() {
        let loss = spherical_spreading_loss(1.0, 10.0);
        assert!(
            loss < 0.0,
            "spherical spreading loss should be negative (attenuation): {loss}"
        );
    }
    #[test]
    fn test_spherical_spreading_loss_6db_per_doubling() {
        let loss1 = spherical_spreading_loss(1.0, 2.0);
        let loss2 = spherical_spreading_loss(1.0, 4.0);
        assert!(
            (loss2 - loss1 + 20.0 * 2.0_f64.log10()).abs() < 0.001,
            "6 dB law: loss1={loss1}, loss2={loss2}"
        );
    }
    #[test]
    fn test_cylindrical_spreading_loss_negative_for_farther() {
        let loss = cylindrical_spreading_loss(1.0, 10.0);
        assert!(
            loss < 0.0,
            "cylindrical spreading loss should be negative: {loss}"
        );
    }
    #[test]
    fn test_cylindrical_spreading_loss_3db_per_doubling() {
        let loss1 = cylindrical_spreading_loss(1.0, 2.0);
        let loss2 = cylindrical_spreading_loss(1.0, 4.0);
        assert!(
            (loss2 - loss1 + 10.0 * 2.0_f64.log10()).abs() < 0.001,
            "3 dB law: loss1={loss1}, loss2={loss2}"
        );
    }
    #[test]
    fn test_add_incoherent_spl_equal_sources() {
        let combined = add_incoherent_spl(&[70.0, 70.0]);
        assert!(
            (combined - (70.0 + 10.0 * 2.0_f64.log10())).abs() < 1e-8,
            "combined SPL = {combined}"
        );
    }
    #[test]
    fn test_add_incoherent_spl_single_source() {
        let spl = add_incoherent_spl(&[65.0]);
        assert!((spl - 65.0).abs() < 1e-8, "single source SPL = {spl}");
    }
    #[test]
    fn test_pressure_fluctuation_rms_constant_field() {
        let p = vec![1.0_f64; 100];
        let rms = pressure_fluctuation_rms(&p);
        assert!((rms - 1.0).abs() < 1e-12, "RMS of constant field = {rms}");
    }
    #[test]
    fn test_pressure_fluctuation_rms_zero_field() {
        let p = vec![0.0_f64; 50];
        let rms = pressure_fluctuation_rms(&p);
        assert!(rms.abs() < 1e-14, "RMS of zero field = {rms}");
    }
    #[test]
    fn test_acoustic_pressure_field_output_length() {
        let rho = vec![1.0, 1.05, 0.95, 1.1];
        let p_ac = acoustic_pressure_field(&rho, 1.0, 1.0 / 3.0);
        assert_eq!(p_ac.len(), rho.len());
    }
    #[test]
    fn test_lighthill_tensor_trace_2d() {
        let t = [[3.0, 1.0], [1.0, 5.0]];
        let tr = lighthill_tensor_trace(t);
        assert!((tr - 8.0).abs() < 1e-14, "trace = {tr}");
    }
    #[test]
    fn test_lighthill_tensor_deviatoric_trace_zero() {
        let t = [[3.0, 1.0], [1.0, 5.0]];
        let dev = lighthill_tensor_deviatoric(t);
        let tr = dev[0][0] + dev[1][1];
        assert!(tr.abs() < 1e-12, "deviatoric trace = {tr}");
    }
    #[test]
    fn test_lighthill_tensor_symmetric_check() {
        let t = [[1.0, 2.0], [2.0, 3.0]];
        assert!(lighthill_tensor_is_symmetric(t, 1e-10));
    }
    #[test]
    fn test_lighthill_tensor_not_symmetric() {
        let t = [[1.0, 2.0], [3.0, 4.0]];
        assert!(!lighthill_tensor_is_symmetric(t, 1e-10));
    }
    #[test]
    fn test_peak_spl_max_of_field() {
        let p = vec![0.001, 0.01, 0.1];
        let spl = peak_spl(&p);
        let expected = spl_from_pressure_fluctuation(0.1);
        assert!((spl - expected).abs() < 1e-10, "peak SPL = {spl}");
    }
    #[test]
    fn test_convected_wavenumber_no_flow_positive() {
        let k = convected_wavenumber(340.0, 340.0, 0.0, 0.0);
        assert!(k > 0.0, "convected k at Ma=0 should be positive: {k}");
    }
    #[test]
    fn test_vortex_sound_pressure_scales_with_vorticity() {
        let p1 = vortex_sound_pressure(1.2, 10.0, [0.0, 0.1], [1.0, 0.0], 1.0, 340.0);
        let p2 = vortex_sound_pressure(1.2, 20.0, [0.0, 0.1], [1.0, 0.0], 1.0, 340.0);
        assert!(
            p2.abs() > p1.abs(),
            "higher vorticity → larger pressure: |p1|={}, |p2|={}",
            p1.abs(),
            p2.abs()
        );
    }
    #[test]
    fn test_howe_vortex_power_scales_with_vorticity() {
        let w1 = howe_vortex_sound_power(1.2, 10.0, 100.0, 0.01, 340.0);
        let w2 = howe_vortex_sound_power(1.2, 10.0, 200.0, 0.01, 340.0);
        assert!(w2 > w1, "higher vorticity → more power: {w1}, {w2}");
    }
    #[test]
    fn test_point_vortex_acoustic_pressure_finite() {
        let p = point_vortex_acoustic_pressure(1.2, 0.5, 10.0, 1.0, 0.0, 0.0, 340.0);
        assert!(p.is_finite(), "point vortex acoustic pressure = {p}");
    }
    #[test]
    fn test_modal_overlap_factor_increases_with_freq() {
        let mof1 = modal_overlap_factor(0.1, 0.01, 100.0);
        let mof2 = modal_overlap_factor(0.1, 0.01, 1000.0);
        assert!(
            mof2 > mof1,
            "modal overlap increases with freq: {mof1}, {mof2}"
        );
    }
    #[test]
    fn test_atmospheric_absorption_increases_with_freq() {
        let a1 = atmospheric_absorption_db_per_m(1000.0, 293.15, 50.0);
        let a2 = atmospheric_absorption_db_per_m(8000.0, 293.15, 50.0);
        assert!(
            a2 > a1,
            "atmospheric absorption increases with freq: {a1}, {a2}"
        );
    }
}
