//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{FwhPanel, FwhSurfacePanel};

/// Compute the Lighthill stress tensor T_ij = rho*u_i*u_j + (p - cs2*rho)*delta_ij.
pub fn lighthill_stress_tensor(rho: f64, u: [f64; 2], p: f64, cs2: f64) -> [[f64; 2]; 2] {
    let excess_p = p - cs2 * rho;
    [
        [rho * u[0] * u[0] + excess_p, rho * u[0] * u[1]],
        [rho * u[1] * u[0], rho * u[1] * u[1] + excess_p],
    ]
}
/// Monopole acoustic source power: Q / (4*pi*rho0*c0).
///
/// `dpdt` is the time derivative of pressure fluctuation (source term Q).
pub fn monopole_source(dpdt: f64, rho0: f64, c0: f64) -> f64 {
    dpdt / (4.0 * PI * rho0 * c0)
}
/// Dipole source strength: |F| / (rho0 * c0).
pub fn dipole_source_strength(force: [f64; 2], rho0: f64, c0: f64) -> f64 {
    let mag = (force[0] * force[0] + force[1] * force[1]).sqrt();
    mag / (rho0 * c0)
}
/// Quadrupole source strength: Frobenius norm of T_ij.
pub fn quadrupole_source_strength(t: [[f64; 2]; 2]) -> f64 {
    let mut sum = 0.0;
    for row in &t {
        for &val in row {
            sum += val * val;
        }
    }
    sum.sqrt()
}
/// Monopole far-field pressure amplitude: p = rho0*c0*Q*k / (4*pi*r) * cos(kr).
///
/// Returns the amplitude of the pressure fluctuation at distance `r`.
pub fn monopole_far_field(q: f64, freq: f64, r: f64, c0: f64, rho0: f64) -> f64 {
    let k = 2.0 * PI * freq / c0;
    rho0 * c0 * q * k / (4.0 * PI * r) * (k * r).cos()
}
/// Dipole far-field pressure: directivity pattern based on observer direction.
///
/// `obs_dir` is the unit vector from source to observer.
pub fn dipole_far_field(
    f: [f64; 2],
    obs_dir: [f64; 2],
    r: f64,
    freq: f64,
    c0: f64,
    rho0: f64,
) -> f64 {
    let k = 2.0 * PI * freq / c0;
    let f_dot_r = f[0] * obs_dir[0] + f[1] * obs_dir[1];
    rho0 * k * f_dot_r / (4.0 * PI * r) * (k * r).cos()
}
/// Sound pressure level in dB: SPL = 20*log10(p_rms / p_ref) where p_ref = 20 µPa.
pub fn sound_pressure_level(p_rms: f64) -> f64 {
    pub(super) const P_REF: f64 = 20e-6;
    20.0 * (p_rms / P_REF).log10()
}
/// Extract acoustic pressure from hydrodynamic pressure.
///
/// p_ac = p_total - c0² * (rho - rho0) for each cell,
/// where rho0 is the mean density.
pub fn acoustic_pressure_from_hydro(p_total: &[f64], rho: &[f64], c0: f64) -> Vec<f64> {
    assert_eq!(
        p_total.len(),
        rho.len(),
        "p_total and rho must have equal length"
    );
    let n = rho.len();
    let rho0 = rho.iter().sum::<f64>() / n as f64;
    let c0sq = c0 * c0;
    p_total
        .iter()
        .zip(rho.iter())
        .map(|(&p, &r)| p - c0sq * (r - rho0))
        .collect()
}
/// Apply Gaussian smoothing to a 2D field stored in row-major order (nx × ny).
///
/// `sigma` controls the smoothing width (in lattice units).
pub fn smooth_acoustic_field(field: &[f64], nx: usize, ny: usize, sigma: f64) -> Vec<f64> {
    assert_eq!(field.len(), nx * ny, "field length must equal nx*ny");
    let radius = (3.0 * sigma).ceil() as isize;
    let mut out = vec![0.0_f64; nx * ny];
    let two_sigma_sq = 2.0 * sigma * sigma;
    for iy in 0..ny {
        for ix in 0..nx {
            let mut weight_sum = 0.0_f64;
            let mut val_sum = 0.0_f64;
            for dy in -radius..=radius {
                for dx in -radius..=radius {
                    let nx_i = ix as isize + dx;
                    let ny_i = iy as isize + dy;
                    if nx_i >= 0 && nx_i < nx as isize && ny_i >= 0 && ny_i < ny as isize {
                        let w = (-(dx * dx + dy * dy) as f64 / two_sigma_sq).exp();
                        weight_sum += w;
                        val_sum += w * field[ny_i as usize * nx + nx_i as usize];
                    }
                }
            }
            out[iy * nx + ix] = val_sum / weight_sum;
        }
    }
    out
}
/// Equilibrium distribution for the acoustic component.
///
/// f_eq = w_i * rho * (1 + (c_i . u) / cs2)  (linearized, acoustic order)
pub fn f_eq_acoustic(rho: f64, u: [f64; 2], w_i: f64, c_i: [f64; 2], cs2: f64) -> f64 {
    let cu = c_i[0] * u[0] + c_i[1] * u[1];
    w_i * rho * (1.0 + cu / cs2)
}
/// Acoustic energy density: E = p_ac²/(2*rho0*c0²) + rho0*|u_ac|²/2.
pub fn acoustic_energy_density(p_ac: f64, u_ac: [f64; 2], rho0: f64, c0: f64) -> f64 {
    let c0sq = c0 * c0;
    let kinetic = rho0 * (u_ac[0] * u_ac[0] + u_ac[1] * u_ac[1]) / 2.0;
    let potential = p_ac * p_ac / (2.0 * rho0 * c0sq);
    potential + kinetic
}
/// Mach number: Ma = velocity / c0.
pub fn mach_number(velocity: f64, c0: f64) -> f64 {
    velocity / c0
}
/// Acoustic wavelength: λ = c0 / freq.
pub fn acoustic_wavelength(freq: f64, c0: f64) -> f64 {
    c0 / freq
}
/// Strouhal number: St = f*L/U.
pub fn strouhal_number(freq: f64, length: f64, velocity: f64) -> f64 {
    freq * length / velocity
}
/// Compute the 3D Lighthill stress tensor T_ij = rho*u_i*u_j + (p - cs2*rho)*delta_ij.
pub fn lighthill_stress_tensor_3d(rho: f64, u: [f64; 3], p: f64, cs2: f64) -> [[f64; 3]; 3] {
    let excess_p = p - cs2 * rho;
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = rho * u[i] * u[j];
            if i == j {
                t[i][j] += excess_p;
            }
        }
    }
    t
}
/// Frobenius norm of a 3x3 tensor: sqrt(sum T_ij^2).
pub fn tensor_frobenius_norm_3d(t: [[f64; 3]; 3]) -> f64 {
    let mut sum = 0.0f64;
    for row in &t {
        for &val in row {
            sum += val * val;
        }
    }
    sum.sqrt()
}
/// Compute the acoustic pressure at an observer from FW-H surface panels.
///
/// Uses the compact-source approximation: loading noise only (no thickness noise).
///
/// `p_ac(x_obs) ~ 1/(4*pi) * sum_panels [ (p_n * n_i * r_hat_i) / r^2 * dA ]`
///
/// where `r` is the distance from panel to observer and `r_hat` is the unit vector.
pub fn fwh_acoustic_pressure(panels: &[FwhPanel], observer: [f64; 3]) -> f64 {
    let mut p_ac = 0.0_f64;
    for panel in panels {
        let rx = observer[0] - panel.center[0];
        let ry = observer[1] - panel.center[1];
        let rz = observer[2] - panel.center[2];
        let r = (rx * rx + ry * ry + rz * rz).sqrt();
        if r < 1e-15 {
            continue;
        }
        let r_hat = [rx / r, ry / r, rz / r];
        let n_dot_r =
            panel.normal[0] * r_hat[0] + panel.normal[1] * r_hat[1] + panel.normal[2] * r_hat[2];
        p_ac += panel.pressure * n_dot_r * panel.area / (4.0 * PI * r * r);
    }
    p_ac
}
/// Monopole directivity: uniform in all directions. Returns 1.0.
pub fn directivity_monopole(_theta: f64) -> f64 {
    1.0
}
/// Dipole directivity: cos(theta).
pub fn directivity_dipole(theta: f64) -> f64 {
    theta.cos()
}
/// Quadrupole directivity: cos²(theta).
pub fn directivity_quadrupole(theta: f64) -> f64 {
    let c = theta.cos();
    c * c
}
/// Lateral quadrupole directivity: sin(theta)*cos(theta).
pub fn directivity_lateral_quadrupole(theta: f64) -> f64 {
    theta.sin() * theta.cos()
}
/// Estimate acoustic power from surface pressure data (monopole radiation).
///
/// `W = sum_panels p_rms^2 * dA / (rho0 * c0)`
pub fn acoustic_power_monopole(p_rms_panels: &[f64], areas: &[f64], rho0: f64, c0: f64) -> f64 {
    assert_eq!(p_rms_panels.len(), areas.len());
    let mut power = 0.0_f64;
    for (prms, area) in p_rms_panels.iter().zip(areas.iter()) {
        power += prms * prms * area / (rho0 * c0);
    }
    power
}
/// Sound power level: PWL = 10*log10(W / W_ref) where W_ref = 10^-12 W.
pub fn sound_power_level(w: f64) -> f64 {
    pub(super) const W_REF: f64 = 1e-12;
    10.0 * (w / W_REF).log10()
}
/// Overall sound pressure level from multiple frequency band SPLs.
///
/// OASPL = 10*log10(sum 10^(SPL_i/10)).
pub fn overall_spl(spl_bands: &[f64]) -> f64 {
    let sum: f64 = spl_bands.iter().map(|&s| 10.0_f64.powf(s / 10.0)).sum();
    10.0 * sum.log10()
}
/// Acoustic intensity: I = p_rms^2 / (rho0 * c0).
pub fn acoustic_intensity(p_rms: f64, rho0: f64, c0: f64) -> f64 {
    p_rms * p_rms / (rho0 * c0)
}
/// Helmholtz number: He = k * L = 2*pi*f*L / c0.
pub fn helmholtz_number(freq: f64, length: f64, c0: f64) -> f64 {
    2.0 * PI * freq * length / c0
}
/// Acoustic impedance: Z = rho0 * c0.
pub fn acoustic_impedance(rho0: f64, c0: f64) -> f64 {
    rho0 * c0
}
/// Near-field pressure amplitude at distance `r` from a monopole.
///
/// Near-field term scales as `1/r^2`:  `p_nf = q * k / (4*pi*r^2)`.
pub fn near_field_pressure_amplitude(q: f64, freq: f64, r: f64, c0: f64) -> f64 {
    let k = 2.0 * PI * freq / c0;
    q * k / (4.0 * PI * r * r)
}
/// Far-field pressure amplitude at distance `r` from a monopole.
///
/// Far-field term scales as `1/r`:  `p_ff = q * k^2 / (4*pi*r)`.
pub fn far_field_pressure_amplitude(q: f64, freq: f64, r: f64, c0: f64) -> f64 {
    let k = 2.0 * PI * freq / c0;
    q * k * k / (4.0 * PI * r)
}
/// Compute the hydrodynamic (incompressible) pressure contribution.
///
/// `p_hydro = -rho * (du/dt * x_i) / r^2`  (simplified estimate from dipole)
pub fn hydrodynamic_pressure(rho: f64, dudt: f64, distance: f64) -> f64 {
    if distance < 1e-15 {
        return 0.0;
    }
    rho * dudt / (distance * distance)
}
/// Extract the acoustic component by subtracting the hydrodynamic pressure.
///
/// `p_ac = p_total - p_hydro`
pub fn acoustic_component_extraction(p_total: f64, p_hydro: f64) -> f64 {
    p_total - p_hydro
}
/// Curle's surface noise analogy.
///
/// Accounts for the scattering of quadrupole sources by solid surfaces.
/// Simplified 2D monopole estimate:
///   `p_ac = dp_dt / (4*pi*r) * cos(omega*(t - r/c0))`
///
/// Returns the amplitude of the surface noise contribution.
pub fn curle_surface_noise(dp_dt: f64, freq: f64, c0: f64, r: f64, t: f64) -> f64 {
    if r < 1e-15 {
        return 0.0;
    }
    let omega = 2.0 * PI * freq;
    dp_dt / (4.0 * PI * r) * (omega * (t - r / c0)).cos()
}
/// Ffowcs Williams-Hall trailing-edge noise analogy.
///
/// Models the sound radiated by turbulence scattering at an edge.
/// Simplified scaling:
///   `p_fwh ∝ u^5/2 * sqrt(rho * L) / (c0^3/2 * r)`
///
/// Returns a scaled amplitude.
pub fn fwh_edge_noise(u: f64, length: f64, c0: f64, r: f64, rho: f64) -> f64 {
    if r < 1e-15 || c0 < 1e-15 {
        return 0.0;
    }
    u.powf(2.5) * (rho * length).sqrt() / (c0.powf(1.5) * r)
}
/// Phillips' quadrupole source decomposition.
///
/// Separates Lighthill's stress tensor into compressibility and entropy terms.
/// Returns the compressibility part: `T_ij_comp = rho * u_i * u_j`.
pub fn phillips_compressibility_term(rho: f64, u_i: f64, u_j: f64) -> f64 {
    rho * u_i * u_j
}
/// Ribner's dilatation noise theory.
///
/// Estimates the acoustic power from isotropic turbulence:
///   `W ~ rho * u^5 / c0^3 * (integral scale)^2`
pub fn ribner_dilatation_noise(rho: f64, u_rms: f64, c0: f64, length_scale: f64) -> f64 {
    rho * u_rms.powi(5) / c0.powi(3) * length_scale * length_scale
}
/// Turbulent boundary layer trailing-edge (TBL-TE) noise.
///
/// Simplified estimate based on Brooks, Pope & Marcolini (BPM) scaling:
///   `p_rms ~ rho * u^{5/2} * delta * / (r * c0^{3/2})`
///
/// Returns the RMS pressure at distance `r`.
pub fn tbl_trailing_edge_noise(
    u: f64,
    rho: f64,
    boundary_layer_thickness: f64,
    r: f64,
    c0: f64,
) -> f64 {
    if r < 1e-15 || c0 < 1e-15 {
        return 0.0;
    }
    rho * u.powf(2.5) * boundary_layer_thickness / (r * c0.powf(1.5))
}
/// Broadband noise scaling with Mach number (5th-power Mach scaling).
///
/// `p ∝ rho * c0^2 * Ma^{5/2}`
pub fn broadband_noise_scaling(mach: f64, rho: f64, c0: f64) -> f64 {
    rho * c0 * c0 * mach.powf(2.5)
}
/// Aeolian tone frequency from Strouhal number (vortex shedding).
///
/// `f_ae = St * U / D` with St ≈ 0.198 (Roshko's value for cylinder).
pub fn aeolian_tone_frequency(velocity: f64, diameter: f64) -> f64 {
    pub(super) const ST: f64 = 0.198;
    ST * velocity / diameter
}
/// Aeolian Strouhal number constant.
pub fn aeolian_strouhal_number() -> f64 {
    0.198
}
/// Jet noise acoustic power (Lighthill's 8th-power law).
///
/// `W = rho * u_j^8 / (c0^5) * A_j * K`
///
/// where `a_j` is the jet exit area and `K` is an empirical constant.
pub fn jet_noise_acoustic_power(u_j: f64, rho: f64, c0: f64, a_j: f64) -> f64 {
    pub(super) const K: f64 = 1e-4;
    rho * u_j.powi(8) / c0.powi(5) * a_j * K
}
/// Coanda-effect noise: interaction of jet with adjacent surface.
///
/// Returns an additional noise term due to surface pressure fluctuations.
pub fn coanda_noise_estimate(u_j: f64, rho: f64, c0: f64, span: f64) -> f64 {
    let ma = u_j / c0;
    rho * c0 * c0 * ma.powi(6) * span
}
/// Turbulent mixing noise level estimate.
///
/// Combines shear layer mixing noise with Lighthill-type quadrupole radiation.
/// Returns dimensionless noise parameter.
pub fn turbulent_mixing_noise_level(u1: f64, u2: f64, rho: f64, c0: f64) -> f64 {
    let du = (u1 - u2).abs();
    let ma = du / c0;
    rho * c0 * c0 * ma * ma * ma * ma * ma
}
/// Acoustic radiation efficiency: ratio of acoustic power to mechanical power.
///
/// `eta = W_ac / (0.5 * rho * u^3 * A)`
pub fn radiation_efficiency(u: f64, c0: f64, rho: f64, area: f64) -> f64 {
    let ma = u / c0;
    let mechanical_power = 0.5 * rho * u.powi(3) * area;
    if mechanical_power < 1e-30 {
        return 0.0;
    }
    let w_ac = rho * c0 * c0 * ma * ma * area;
    w_ac / mechanical_power
}
/// Extract plane-wave mode amplitude from a ring array of pressures.
///
/// For a uniform incident plane wave all sensors give the same value,
/// so the plane-wave mode = mean of the pressures.
pub fn plane_wave_mode(pressures: &[f64]) -> f64 {
    if pressures.is_empty() {
        return 0.0;
    }
    pressures.iter().sum::<f64>() / pressures.len() as f64
}
/// Atmospheric sound attenuation coefficient (ISO 9613-1 simplified).
///
/// Returns the attenuation coefficient `alpha` \[dB/m\] as a function of
/// frequency `f` \[Hz\], relative humidity `hr` \[%\], and temperature `T` \[°C\].
///
/// Simplified formula (valid for moderate conditions):
///   `alpha ≈ 8.686 * f^2 * (1.84e-11 / (hr * T_K^0.5) + 0.1exp(-T_K/2239.1))`
pub fn atmospheric_attenuation(freq: f64, humidity_pct: f64, temp_c: f64) -> f64 {
    let t_k = temp_c + 273.15;
    let term1 = 1.84e-11 * freq * freq / (humidity_pct * t_k.sqrt());
    let term2 = 0.1 * (-t_k / 2239.1_f64).exp() * freq * freq / t_k;
    8.686 * (term1 + term2)
}
/// Reflection coefficient at a hard (rigid) wall.
///
/// For a rigid surface, `R = +1` (pressure doubles at the wall).
/// Returns the incident pressure (equal contribution from image source).
pub fn hard_wall_reflection(p_incident: f64) -> f64 {
    p_incident
}
/// Sound absorption coefficient from reflection coefficient.
///
/// `alpha_abs = 1 - |R|^2`
///
/// `r_real`, `r_imag` are the real and imaginary parts of the reflection
/// coefficient.
pub fn sound_absorption_coefficient(r_real: f64, r_imag: f64) -> f64 {
    let r_mag_sq = r_real * r_real + r_imag * r_imag;
    (1.0 - r_mag_sq).clamp(0.0, 1.0)
}
/// Maekawa barrier attenuation (simplified).
///
/// `A_bar ≈ 10 * log10(3 + 20 * N)` where `N = 2*delta / lambda`
/// is the Fresnel number, `delta` is the path length difference \[m\],
/// and `lambda` is the acoustic wavelength.
pub fn maekawa_barrier_attenuation(delta_path: f64, r: f64, c0: f64) -> f64 {
    let freq = c0 / r;
    let lambda = c0 / freq;
    let fresnel_n = 2.0 * delta_path / lambda;
    10.0 * (3.0 + 20.0 * fresnel_n).log10()
}
/// Image source pressure at a reflecting wall.
///
/// Adds the contribution of the mirror image source.  For a hard wall at
/// the midpoint, observer and image are equidistant:
///   `p_total = p_direct + p_image = 2 * p(r)`
///
/// Here we use the far-field monopole formula for both contributions.
pub fn image_source_pressure(q: f64, freq: f64, c0: f64, r: f64) -> f64 {
    let p = monopole_far_field(q, freq, r, c0, 1.0);
    p + p
}
/// Generalized Lighthill source term including entropy fluctuations.
///
/// `T_ij = rho * u_i * u_j + (p - cs2*rho)*delta_ij - tau_ij`
///
/// where `tau_ij` is the viscous stress tensor.  This function returns the
/// isotropic part `p - cs2 * rho - tau_ii / 3`.
pub fn lighthill_entropy_term(p: f64, cs2: f64, rho: f64, tau_trace: f64) -> f64 {
    p - cs2 * rho - tau_trace / 3.0
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::aeroacoustics::types::*;
    #[test]
    fn test_lighthill_diagonal_excess_pressure() {
        let rho = 1.2;
        let u = [0.1, 0.0];
        let p = 0.4;
        let cs2 = 1.0 / 3.0;
        let t = lighthill_stress_tensor(rho, u, p, cs2);
        let expected_00 = rho * u[0] * u[0] + (p - cs2 * rho);
        assert!((t[0][0] - expected_00).abs() < 1e-14);
    }
    #[test]
    fn test_lighthill_symmetry() {
        let t = lighthill_stress_tensor(1.0, [0.05, 0.03], 0.35, 1.0 / 3.0);
        assert!((t[0][1] - t[1][0]).abs() < 1e-14);
    }
    #[test]
    fn test_lighthill_zero_velocity() {
        let rho = 1.0;
        let p = 1.0 / 3.0;
        let cs2 = 1.0 / 3.0;
        let t = lighthill_stress_tensor(rho, [0.0, 0.0], p, cs2);
        assert!(t[0][1].abs() < 1e-14);
        assert!(t[1][0].abs() < 1e-14);
        assert!(t[0][0].abs() < 1e-14);
        assert!(t[1][1].abs() < 1e-14);
    }
    #[test]
    fn test_monopole_source_scaling() {
        let rho0 = 1.2;
        let c0 = 340.0;
        let s1 = monopole_source(1.0, rho0, c0);
        let s2 = monopole_source(2.0, rho0, c0);
        assert!((s2 / s1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_monopole_source_zero() {
        assert!(monopole_source(0.0, 1.2, 340.0).abs() < 1e-14);
    }
    #[test]
    fn test_dipole_source_strength() {
        let rho0 = 1.2;
        let c0 = 340.0;
        let d = dipole_source_strength([3.0, 4.0], rho0, c0);
        let expected = 5.0 / (rho0 * c0);
        assert!((d - expected).abs() < 1e-12);
    }
    #[test]
    fn test_dipole_zero_force() {
        assert!(dipole_source_strength([0.0, 0.0], 1.0, 340.0).abs() < 1e-14);
    }
    #[test]
    fn test_quadrupole_frobenius_norm() {
        let t = [[1.0, 0.0], [0.0, 1.0]];
        let q = quadrupole_source_strength(t);
        assert!((q - 2.0_f64.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_quadrupole_zero() {
        assert!(quadrupole_source_strength([[0.0, 0.0], [0.0, 0.0]]).abs() < 1e-14);
    }
    #[test]
    fn test_monopole_far_field_value() {
        let q = 1.0;
        let freq = 100.0;
        let c0 = 340.0;
        let rho0 = 1.2;
        let r = 10.0;
        let k = 2.0 * PI * freq / c0;
        let expected = rho0 * c0 * q * k / (4.0 * PI * r) * (k * r).cos();
        let result = monopole_far_field(q, freq, r, c0, rho0);
        assert!(
            (result - expected).abs() < 1e-12,
            "Monopole far-field mismatch"
        );
    }
    #[test]
    fn test_spl_reference_pressure() {
        let spl = sound_pressure_level(20e-6);
        assert!(
            spl.abs() < 1e-10,
            "SPL at reference should be 0 dB, got {spl}"
        );
    }
    #[test]
    fn test_spl_ten_times_pressure() {
        let spl = sound_pressure_level(200e-6);
        assert!((spl - 20.0).abs() < 1e-10, "SPL = {spl}, expected 20 dB");
    }
    #[test]
    fn test_acoustic_pressure_uniform() {
        let rho0 = 1.0_f64;
        let c0 = (1.0_f64 / 3.0_f64).sqrt();
        let c0sq = c0 * c0;
        let p = vec![c0sq * rho0; 9];
        let rho = vec![rho0; 9];
        let p_ac = acoustic_pressure_from_hydro(&p, &rho, c0);
        for &val in &p_ac {
            assert!((val - c0sq * rho0).abs() < 1e-12, "Uniform: p_ac = {val}");
        }
    }
    #[test]
    fn test_acoustic_pressure_non_uniform() {
        let rho = vec![1.0, 1.0, 1.0, 1.0, 1.5];
        let c0 = 1.0;
        let c0sq = c0 * c0;
        let rho0 = rho.iter().sum::<f64>() / rho.len() as f64;
        let p: Vec<f64> = rho.iter().map(|&r| c0sq * r).collect();
        let p_ac = acoustic_pressure_from_hydro(&p, &rho, c0);
        for &val in &p_ac {
            assert!((val - c0sq * rho0).abs() < 1e-12, "p_ac = {val}");
        }
    }
    #[test]
    fn test_smooth_uniform_field() {
        let nx = 5;
        let ny = 5;
        let field = vec![2.5_f64; nx * ny];
        let smoothed = smooth_acoustic_field(&field, nx, ny, 1.0);
        for &v in &smoothed {
            assert!((v - 2.5).abs() < 1e-10, "Smoothed uniform = {v}");
        }
    }
    #[test]
    fn test_smooth_spike() {
        let nx = 7;
        let ny = 7;
        let mut field = vec![0.0_f64; nx * ny];
        field[3 * nx + 3] = 1.0;
        let smoothed = smooth_acoustic_field(&field, nx, ny, 1.0);
        assert!(
            smoothed[3 * nx + 4] > 0.0,
            "Neighbor should be non-zero after smoothing"
        );
        assert!(
            smoothed[3 * nx + 3] < 1.0,
            "Center should be reduced after smoothing"
        );
    }
    #[test]
    fn test_f_eq_acoustic_zero_velocity() {
        let rho = 1.5;
        let w = 4.0 / 9.0;
        let cs2 = 1.0 / 3.0;
        let f = f_eq_acoustic(rho, [0.0, 0.0], w, [0.0, 0.0], cs2);
        assert!((f - w * rho).abs() < 1e-14);
    }
    #[test]
    fn test_acoustic_energy_potential_only() {
        let p_ac = 1.0;
        let rho0 = 1.2;
        let c0 = 340.0;
        let e = acoustic_energy_density(p_ac, [0.0, 0.0], rho0, c0);
        let expected = p_ac * p_ac / (2.0 * rho0 * c0 * c0);
        assert!((e - expected).abs() < 1e-12);
    }
    #[test]
    fn test_acoustic_energy_kinetic_only() {
        let rho0 = 1.2;
        let u_ac = [1.0, 0.0];
        let e = acoustic_energy_density(0.0, u_ac, rho0, 340.0);
        let expected = rho0 * 0.5;
        assert!((e - expected).abs() < 1e-12);
    }
    #[test]
    fn test_helmholtz_frequency_positive() {
        let hr = HelmholtzResonator::new(0.001, 0.05, 0.0001, 343.0);
        let f = hr.resonance_frequency();
        assert!(f > 0.0, "Resonance frequency should be positive: {f}");
    }
    #[test]
    fn test_helmholtz_frequency_volume_effect() {
        let hr_small = HelmholtzResonator::new(0.0005, 0.05, 0.0001, 343.0);
        let hr_large = HelmholtzResonator::new(0.002, 0.05, 0.0001, 343.0);
        let f_small = hr_small.resonance_frequency();
        let f_large = hr_large.resonance_frequency();
        assert!(
            f_large < f_small,
            "Larger volume should give lower frequency"
        );
    }
    #[test]
    fn test_helmholtz_transmission_loss_positive() {
        let hr = HelmholtzResonator::new(0.001, 0.05, 0.0001, 343.0);
        let tl = hr.transmission_loss_peak();
        assert!(tl > 0.0, "Transmission loss should be positive: {tl}");
    }
    #[test]
    fn test_mach_number_subsonic() {
        let ma = mach_number(100.0, 340.0);
        assert!(ma < 1.0 && ma > 0.0, "Mach = {ma}");
    }
    #[test]
    fn test_mach_number_sonic() {
        let ma = mach_number(340.0, 340.0);
        assert!((ma - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_wavelength_inverse_frequency() {
        let lambda1 = acoustic_wavelength(100.0, 340.0);
        let lambda2 = acoustic_wavelength(200.0, 340.0);
        assert!((lambda2 * 2.0 - lambda1).abs() < 1e-12);
    }
    #[test]
    fn test_strouhal_number() {
        let st = strouhal_number(10.0, 0.1, 5.0);
        assert!((st - 0.2).abs() < 1e-14, "St = {st}, expected 0.2");
    }
    #[test]
    fn test_lighthill_3d_diagonal() {
        let rho = 1.2;
        let u = [0.1, 0.0, 0.0];
        let p = 0.4;
        let cs2 = 1.0 / 3.0;
        let t = lighthill_stress_tensor_3d(rho, u, p, cs2);
        let expected_00 = rho * u[0] * u[0] + (p - cs2 * rho);
        assert!((t[0][0] - expected_00).abs() < 1e-14);
    }
    #[test]
    fn test_lighthill_3d_symmetry() {
        let t = lighthill_stress_tensor_3d(1.0, [0.05, 0.03, 0.01], 0.35, 1.0 / 3.0);
        for (i, row) in t.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - t[j][i]).abs() < 1e-14,
                    "3D Lighthill not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_lighthill_3d_zero_velocity() {
        let t = lighthill_stress_tensor_3d(1.0, [0.0, 0.0, 0.0], 1.0 / 3.0, 1.0 / 3.0);
        for (i, row) in t.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if i != j {
                    assert!(val.abs() < 1e-14);
                }
            }
        }
    }
    #[test]
    fn test_frobenius_norm_3d() {
        let t = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let norm = tensor_frobenius_norm_3d(t);
        assert!((norm - 3.0_f64.sqrt()).abs() < 1e-12, "norm = {norm}");
    }
    #[test]
    fn test_fwh_single_panel() {
        let panel = FwhPanel {
            center: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            area: 1.0,
            pressure: 1.0,
            un: 0.0,
        };
        let observer = [10.0, 0.0, 0.0];
        let p_ac = fwh_acoustic_pressure(&[panel], observer);
        assert!(p_ac > 0.0, "FW-H p_ac should be positive: {p_ac}");
    }
    #[test]
    fn test_fwh_perpendicular_panel() {
        let panel = FwhPanel {
            center: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            area: 1.0,
            pressure: 1.0,
            un: 0.0,
        };
        let observer = [10.0, 0.0, 0.0];
        let p_ac = fwh_acoustic_pressure(&[panel], observer);
        assert!(p_ac.abs() < 1e-14, "FW-H perpendicular: p_ac = {p_ac}");
    }
    #[test]
    fn test_directivity_monopole() {
        assert!((directivity_monopole(0.0) - 1.0).abs() < 1e-14);
        assert!((directivity_monopole(PI / 2.0) - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_directivity_dipole_90() {
        let d = directivity_dipole(PI / 2.0);
        assert!(d.abs() < 1e-14, "dipole at 90°: {d}");
    }
    #[test]
    fn test_directivity_quadrupole_90() {
        let d = directivity_quadrupole(PI / 2.0);
        assert!(d.abs() < 1e-14, "quadrupole at 90°: {d}");
    }
    #[test]
    fn test_acoustic_power_zero() {
        let w = acoustic_power_monopole(&[0.0, 0.0], &[1.0, 1.0], 1.2, 340.0);
        assert!(w.abs() < 1e-20, "power = {w}");
    }
    #[test]
    fn test_acoustic_power_positive() {
        let w = acoustic_power_monopole(&[1.0], &[1.0], 1.2, 340.0);
        assert!(w > 0.0, "power should be positive: {w}");
    }
    #[test]
    fn test_sound_power_level_scaling() {
        let w1 = 1e-6;
        let w2 = 2e-6;
        let pwl1 = sound_power_level(w1);
        let pwl2 = sound_power_level(w2);
        let diff = pwl2 - pwl1;
        assert!(
            (diff - 10.0 * 2.0_f64.log10()).abs() < 1e-10,
            "PWL diff = {diff}"
        );
    }
    #[test]
    fn test_overall_spl_single_band() {
        let oaspl = overall_spl(&[80.0]);
        assert!((oaspl - 80.0).abs() < 1e-10, "OASPL = {oaspl}");
    }
    #[test]
    fn test_overall_spl_two_equal_bands() {
        let oaspl = overall_spl(&[80.0, 80.0]);
        let expected = 80.0 + 10.0 * 2.0_f64.log10();
        assert!(
            (oaspl - expected).abs() < 1e-10,
            "OASPL = {oaspl}, expected {expected}"
        );
    }
    #[test]
    fn test_acoustic_intensity() {
        let i = acoustic_intensity(2.0, 1.2, 340.0);
        let expected = 4.0 / (1.2 * 340.0);
        assert!(
            (i - expected).abs() < 1e-12,
            "intensity = {i}, expected {expected}"
        );
    }
    #[test]
    fn test_helmholtz_number() {
        let he = helmholtz_number(340.0, 1.0, 340.0);
        assert!((he - 2.0 * PI).abs() < 1e-10, "He = {he}");
    }
    #[test]
    fn test_acoustic_impedance() {
        let z = acoustic_impedance(1.2, 340.0);
        assert!((z - 408.0).abs() < 1e-10, "Z = {z}");
    }
    #[test]
    fn test_directivity_lateral_quadrupole() {
        let d = directivity_lateral_quadrupole(PI / 4.0);
        assert!((d - 0.5).abs() < 1e-10, "lateral quad at 45°: {d}");
    }
    #[test]
    fn test_near_field_dominates_at_short_range() {
        let nf = near_field_pressure_amplitude(1.0, 1.0, 10.0, 340.0);
        let ff = far_field_pressure_amplitude(1.0, 1.0, 10.0, 340.0);
        assert!(
            nf > ff * 0.9,
            "near-field should dominate at r=1: nf={nf}, ff={ff}"
        );
    }
    #[test]
    fn test_far_field_dominates_at_large_range() {
        let nf = near_field_pressure_amplitude(1.0, 1.0, 100.0, 340.0);
        let ff = far_field_pressure_amplitude(1.0, 1.0, 100.0, 340.0);
        assert!(
            ff >= nf * 0.9,
            "far-field should dominate at r=100: ff={ff}, nf={nf}"
        );
    }
    #[test]
    fn test_near_field_inverse_square_decay() {
        let p1 = near_field_pressure_amplitude(1.0, 1.0, 5.0, 340.0);
        let p2 = near_field_pressure_amplitude(1.0, 1.0, 10.0, 340.0);
        let ratio = p1 / p2;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "NF: 2x distance → 4x decay, ratio={ratio}"
        );
    }
    #[test]
    fn test_far_field_inverse_decay() {
        let p1 = far_field_pressure_amplitude(1.0, 1.0, 5.0, 340.0);
        let p2 = far_field_pressure_amplitude(1.0, 1.0, 10.0, 340.0);
        let ratio = p1 / p2;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "FF: 2x distance → 2x decay, ratio={ratio}"
        );
    }
    #[test]
    fn test_curle_surface_noise_positive() {
        let p_surf = curle_surface_noise(1.0, 1.0, 340.0, 10.0, 0.0);
        assert!(
            p_surf.is_finite(),
            "Curle surface noise should be finite: {p_surf}"
        );
    }
    #[test]
    fn test_curle_zero_surface_pressure() {
        let p_surf = curle_surface_noise(0.0, 1.0, 340.0, 10.0, 0.0);
        assert!(p_surf.abs() < 1e-20, "Zero surface pressure: {p_surf}");
    }
    #[test]
    fn test_fwh_edge_noise_positive() {
        let p = fwh_edge_noise(1.0, 0.1, 340.0, 10.0, 1.0);
        assert!(p >= 0.0, "FWH edge noise should be non-negative: {p}");
    }
    #[test]
    fn test_tbl_te_noise_positive() {
        let spl = tbl_trailing_edge_noise(10.0, 1.0, 0.01, 1.0, 340.0);
        assert!(spl > 0.0, "TBL-TE noise SPL should be positive: {spl}");
    }
    #[test]
    fn test_broadband_noise_mach_scaling() {
        let p1 = broadband_noise_scaling(0.1, 1.2, 340.0);
        let p2 = broadband_noise_scaling(0.2, 1.2, 340.0);
        assert!(
            p2 > p1,
            "Higher Mach → more broadband noise: p1={p1}, p2={p2}"
        );
    }
    #[test]
    fn test_aeolian_tone_frequency() {
        let f = aeolian_tone_frequency(10.0, 0.01);
        let expected = 0.198 * 10.0 / 0.01;
        assert!(
            (f - expected).abs() < 1e-10,
            "Aeolian frequency = {f}, expected {expected}"
        );
    }
    #[test]
    fn test_aeolian_strouhal_number() {
        let st = aeolian_strouhal_number();
        assert!((st - 0.198).abs() < 1e-10, "Aeolian St = {st}");
    }
    #[test]
    fn test_jet_noise_eighth_power() {
        let w1 = jet_noise_acoustic_power(10.0, 1.2, 340.0, 1e-4);
        let w2 = jet_noise_acoustic_power(20.0, 1.2, 340.0, 1e-4);
        let ratio = w2 / w1;
        assert!(
            (ratio - 256.0).abs() / 256.0 < 0.01,
            "ratio = {ratio}, expected ~256"
        );
    }
    #[test]
    fn test_radiation_efficiency_positive() {
        let eta = radiation_efficiency(1.0, 340.0, 1.2, 1e-2);
        assert!(
            eta >= 0.0,
            "Radiation efficiency should be non-negative: {eta}"
        );
    }
    #[test]
    fn test_plane_wave_mode_extraction() {
        let pressures = vec![1.0f64; 4];
        let p_mode = plane_wave_mode(&pressures);
        assert!(
            (p_mode - 1.0).abs() < 1e-12,
            "Uniform field → plane wave = 1.0: {p_mode}"
        );
    }
    #[test]
    fn test_plane_wave_mode_cancellation() {
        let pressures = vec![1.0, -1.0, 1.0, -1.0];
        let p_mode = plane_wave_mode(&pressures);
        assert!(p_mode.abs() < 1e-12, "Alternating: plane wave = {p_mode}");
    }
    #[test]
    fn test_attenuation_coefficient_positive() {
        let alpha = atmospheric_attenuation(1000.0, 50.0, 20.0);
        assert!(alpha > 0.0, "Attenuation should be positive: {alpha}");
    }
    #[test]
    fn test_attenuation_increases_with_freq() {
        let a1 = atmospheric_attenuation(1000.0, 50.0, 20.0);
        let a2 = atmospheric_attenuation(4000.0, 50.0, 20.0);
        assert!(
            a2 > a1,
            "Higher frequency → more attenuation: a1={a1}, a2={a2}"
        );
    }
    #[test]
    fn test_hard_wall_reflection() {
        let p_inc = 1.0;
        let p_refl = hard_wall_reflection(p_inc);
        assert!(
            (p_refl - 1.0).abs() < 1e-12,
            "Hard wall reflection = {p_refl}"
        );
    }
    #[test]
    fn test_sound_absorption_coefficient_range() {
        let alpha = sound_absorption_coefficient(0.5, 0.5);
        assert!(
            (0.0..=1.0).contains(&alpha),
            "absorption should be in [0,1]: {alpha}"
        );
    }
    #[test]
    fn test_maekawa_attenuation_positive() {
        let att = maekawa_barrier_attenuation(1.0, 10.0, 340.0);
        assert!(
            att >= 0.0,
            "Maekawa attenuation should be non-negative: {att}"
        );
    }
    #[test]
    fn test_image_source_doubling() {
        let q = 1.0;
        let freq = 1.0;
        let c0 = 340.0;
        let r = 1.0;
        let p_direct = monopole_far_field(q, freq, r, c0, 1.0);
        let p_total = image_source_pressure(q, freq, c0, r);
        assert!(
            (p_total - 2.0 * p_direct).abs() < 1e-12,
            "Image source should double: p_total={p_total}, 2*p_direct={}",
            2.0 * p_direct
        );
    }
}
/// Compute FWH thickness noise contribution from one panel.
///
/// p'_T = ρ_0 * (v_n * dA) / (4π r c_0) * r̂ · n̂
pub fn fwh_thickness_noise(panel: &FwhSurfacePanel, rho0: f64, c0: f64) -> f64 {
    let v_n = dot3(panel.velocity, panel.normal);
    let r_dot_n = dot3(panel.r_hat, panel.normal);
    rho0 * v_n * panel.area * r_dot_n / (4.0 * PI * panel.r * c0)
}
/// Compute FWH loading noise contribution from one panel.
///
/// p'_L = Δp * (r̂ · n̂) * dA / (4π r c_0)
pub fn fwh_loading_noise(panel: &FwhSurfacePanel, c0: f64) -> f64 {
    let r_dot_n = dot3(panel.r_hat, panel.normal);
    panel.delta_p * r_dot_n * panel.area / (4.0 * PI * panel.r * c0)
}
/// Sum FWH contributions over all panels to get total acoustic pressure.
pub fn fwh_total_pressure(panels: &[FwhSurfacePanel], rho0: f64, c0: f64) -> f64 {
    panels
        .iter()
        .map(|p| fwh_thickness_noise(p, rho0, c0) + fwh_loading_noise(p, c0))
        .sum()
}
/// Near-field acoustic pressure from a monopole source.
///
/// Includes both far-field (1/r) and near-field (1/r^2) terms.
/// p(r,t) = Q/(4π) * \[1/(r c_0) * dpdt + 1/r^2 * p_static\]
pub fn near_field_monopole_pressure(q: f64, dpdt: f64, p_static: f64, r: f64, c0: f64) -> f64 {
    if r < 1e-15 {
        return 0.0;
    }
    q / (4.0 * PI) * (dpdt / (r * c0) + p_static / (r * r))
}
/// Dipole near-field pressure along the dipole axis.
///
/// p = F_0 * cos(θ) / (4π) * \[1/(r c_0) * dFdt/F_0 + 1/r^2\]
pub fn near_field_dipole_pressure(
    force_magnitude: f64,
    dfdt: f64,
    r: f64,
    cos_theta: f64,
    c0: f64,
) -> f64 {
    if r < 1e-15 {
        return 0.0;
    }
    force_magnitude * cos_theta / (4.0 * PI)
        * (dfdt / (force_magnitude.max(1e-30) * r * c0) + 1.0 / (r * r))
}
/// Reactive near-field intensity (imaginary part of acoustic intensity).
///
/// I_reactive = p'^2 / (2 ρ_0 ω r) (near-field evanescent contribution)
pub fn reactive_near_field_intensity(p_rms: f64, rho0: f64, omega: f64, r: f64) -> f64 {
    if rho0 < 1e-15 || omega < 1e-15 || r < 1e-15 {
        return 0.0;
    }
    p_rms * p_rms / (2.0 * rho0 * omega * r)
}
/// Broadband trailing-edge noise spectral density (Brooks-Pope-Marcolini model, simplified).
///
/// Returns the sound pressure level contribution in dB at frequency `f`.
/// Scaling: SPL = 10*log10(d * M^5 * r_e^{-2} * H(St))
/// where St = f*d/U is Strouhal number based on displacement thickness `d`.
pub fn trailing_edge_noise_spl(f: f64, u: f64, d: f64, r_e: f64, rho0: f64, c0: f64) -> f64 {
    if u < 1e-15 || d < 1e-15 || r_e < 1e-15 {
        return f64::NEG_INFINITY;
    }
    let m = u / c0;
    let st = f * d / u;
    let h = (-(st - 0.2).powi(2) / (2.0 * 0.05_f64.powi(2))).exp();
    let p_ref = 20e-6_f64;
    let p2 = rho0 * c0.powi(2) * d * m.powi(5) * h / (r_e * r_e);
    if p2 <= 0.0 {
        return f64::NEG_INFINITY;
    }
    10.0 * (p2 / (p_ref * p_ref)).log10()
}
/// A-weighted sound pressure level correction (dB) at frequency f (Hz).
///
/// Uses the standard A-weighting formula.
pub fn a_weighting_db(f: f64) -> f64 {
    if f < 1.0 {
        return -100.0;
    }
    let f2 = f * f;
    let num = 12194.0_f64.powi(2) * f2 * f2;
    let den = (f2 + 20.6_f64.powi(2))
        * ((f2 + 107.7_f64.powi(2)) * (f2 + 737.9_f64.powi(2))).sqrt()
        * (f2 + 12194.0_f64.powi(2));
    if den < 1e-30 {
        return -100.0;
    }
    2.0 + 20.0 * (num / den).log10()
}
/// Compute overall A-weighted SPL from a spectrum.
///
/// `freqs` and `spls` must have the same length (dB values per frequency band).
pub fn overall_a_weighted_spl(freqs: &[f64], spls: &[f64]) -> f64 {
    if freqs.is_empty() {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = freqs
        .iter()
        .zip(spls.iter())
        .map(|(&f, &spl)| {
            let spl_a = spl + a_weighting_db(f);
            10.0_f64.powf(spl_a / 10.0)
        })
        .sum();
    if sum <= 0.0 {
        return f64::NEG_INFINITY;
    }
    10.0 * sum.log10()
}
/// Aeolian tone frequency from a circular cylinder (Strouhal law).
///
/// f_s = St * U / D
/// St ≈ 0.2 for Re in 300..3×10^5.
pub fn aeolian_tone_freq_strouhal(u: f64, d: f64, strouhal: f64) -> f64 {
    strouhal * u / d
}
/// Strouhal number for a cylinder from flow velocity, frequency, and diameter.
pub fn cylinder_strouhal(f_s: f64, u: f64, d: f64) -> f64 {
    if u < 1e-15 || d < 1e-15 {
        return 0.0;
    }
    f_s * d / u
}
/// Lift-fluctuation tonal noise SPL estimate (Curle's analogy, 2D cross-section).
///
/// SPL ≈ 10*log10( F_L^2 / (4π^2 r^2 ρ_0 c_0^2 p_ref^2) * (ω/c_0)^2 )
pub fn tonal_lift_noise_spl(f_lift: f64, omega: f64, r: f64, rho0: f64, c0: f64) -> f64 {
    let p_ref = 20e-6_f64;
    let denom = 4.0 * PI * PI * r * r * rho0 * c0 * c0 * p_ref * p_ref;
    if denom < 1e-30 {
        return f64::NEG_INFINITY;
    }
    let p2 = f_lift * f_lift * omega * omega / (c0 * c0) / denom;
    if p2 <= 0.0 {
        return f64::NEG_INFINITY;
    }
    10.0 * p2.log10()
}
/// Lock-in bandwidth: range of reduced velocity U* = U/(f_s * D) over which
/// vortex shedding locks to body oscillation.
///
/// Returns (U*_lower, U*_upper) based on empirical ±20% of nominal St.
pub fn lock_in_bandwidth(strouhal_nominal: f64) -> (f64, f64) {
    let u_star_nominal = 1.0 / strouhal_nominal;
    (u_star_nominal * 0.80, u_star_nominal * 1.20)
}
/// Time-averaged acoustic intensity (W/m^2): I = p'^2 / (2 ρ_0 c_0).
pub fn acoustic_intensity_rms(p_rms: f64, rho0: f64, c0: f64) -> f64 {
    p_rms * p_rms / (2.0 * rho0 * c0)
}
/// Acoustic intensity vector in direction `r_hat` (unit vector).
///
/// I_vec = I * r_hat
pub fn acoustic_intensity_rms_vector(p_rms: f64, rho0: f64, c0: f64, r_hat: [f64; 3]) -> [f64; 3] {
    let i = acoustic_intensity_rms(p_rms, rho0, c0);
    [i * r_hat[0], i * r_hat[1], i * r_hat[2]]
}
/// Acoustic energy density: E = p'^2 / (ρ_0 c_0^2).
pub fn acoustic_energy_density_rms(p_rms: f64, rho0: f64, c0: f64) -> f64 {
    p_rms * p_rms / (rho0 * c0 * c0)
}
/// Active acoustic power through a surface of area A: W = I * A.
pub fn acoustic_power(p_rms: f64, rho0: f64, c0: f64, area: f64) -> f64 {
    acoustic_intensity_rms(p_rms, rho0, c0) * area
}
/// Sound power level: L_W = 10 * log10(W / W_ref)  \[dB re 1e-12 W\].
pub fn sound_power_level_db(w: f64) -> f64 {
    let w_ref = 1e-12_f64;
    if w <= 0.0 {
        return f64::NEG_INFINITY;
    }
    10.0 * (w / w_ref).log10()
}
/// Convert sound power level (dBW re 1e-12 W) back to watts.
pub fn sound_power_from_level(lw_db: f64) -> f64 {
    1e-12_f64 * 10.0_f64.powf(lw_db / 10.0)
}
/// Directivity index: DI = 10 * log10(I_dir / I_omni).
///
/// `i_directional` is the intensity in the direction of interest,
/// `total_power` is the total radiated power, `r` is the measurement distance.
pub fn directivity_index(i_directional: f64, total_power: f64, r: f64) -> f64 {
    let i_omni = total_power / (4.0 * PI * r * r);
    if i_omni < 1e-30 {
        return 0.0;
    }
    10.0 * (i_directional / i_omni).log10()
}
/// Radiated sound power from a vibrating sphere (simplified).
///
/// W = 0.5 * ρ_0 * c_0 * (k*a)^2 / (1 + (k*a)^2) * v_rms^2 * 4π a^2
/// where k = ω/c_0, a = sphere radius.
pub fn sphere_radiated_power(v_rms: f64, omega: f64, a: f64, rho0: f64, c0: f64) -> f64 {
    let k = omega / c0;
    let ka = k * a;
    let radiation_efficiency = ka * ka / (1.0 + ka * ka);
    let area = 4.0 * PI * a * a;
    0.5 * rho0 * c0 * radiation_efficiency * v_rms * v_rms * area
}
/// Dot product of two 3D vectors.
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// 2D free-space acoustic Green's function magnitude: G = ln(1/(k*r)) / (2π).
///
/// Returns the real part of the 2D monopole Green's function at wavenumber `k`
/// and radial distance `r` (r > 0).
pub fn greens_function_2d(k: f64, r: f64) -> f64 {
    if r <= 0.0 {
        return 0.0;
    }
    -(k * r).ln() / (2.0 * PI)
}
/// 3D free-space acoustic Green's function: G = exp(-i k r) / (4π r).
///
/// Returns the amplitude |G| = 1 / (4π r).
pub fn greens_function_3d_amplitude(r: f64) -> f64 {
    if r <= 0.0 {
        return 0.0;
    }
    1.0 / (4.0 * PI * r)
}
/// Half-space Green's function using a single image source.
///
/// p = G(r_direct) + G(r_image) where the image source is the reflection
/// of the original source across the wall (y = 0 plane).
pub fn half_space_greens(r_direct: f64, r_image: f64) -> f64 {
    greens_function_3d_amplitude(r_direct) + greens_function_3d_amplitude(r_image)
}
/// Compute the retarded time for a moving source: t_ret = t - r / c0.
pub fn retarded_time(t: f64, r: f64, c0: f64) -> f64 {
    t - r / c0
}
/// Plane-wave scattering by a rigid sphere (first-order correction).
///
/// Returns the scattered pressure amplitude at angle `theta` (rad) and
/// distance `r` from sphere of radius `a`:
/// p_sc ≈ -i * (k*a)³ / 3 * cos(theta) * (1/(k*r)) * p_inc
/// (Rayleigh limit: k*a << 1)
pub fn sphere_scattering_rayleigh(p_inc: f64, k: f64, a: f64, r: f64, theta: f64) -> f64 {
    let ka = k * a;
    let ka3 = ka * ka * ka;
    ka3 / 3.0 * theta.cos().abs() / (k * r) * p_inc.abs()
}
/// Cylinder scattering in 2D: monopole scattered amplitude.
///
/// In the Rayleigh limit: p_sc ≈ (π/4) * (ka)² * p_inc / sqrt(kr)
pub fn cylinder_scattering_rayleigh(p_inc: f64, k: f64, a: f64, r: f64) -> f64 {
    let ka = k * a;
    PI / 4.0 * ka * ka * p_inc.abs() / (k * r).sqrt()
}
/// Total scattering cross-section for a rigid sphere in the Rayleigh limit.
///
/// σ_scat = (2π/3) * (k*a)^4 * a²  (valid for k*a << 1)
pub fn sphere_scattering_cross_section(k: f64, a: f64) -> f64 {
    let ka = k * a;
    (2.0 * PI / 3.0) * ka * ka * ka * ka * a * a
}
/// Diffraction by a semi-infinite rigid wedge (Sommerfeld solution, approximation).
///
/// Returns the SPL correction dB relative to free-field.
pub fn wedge_diffraction_spl_correction(delta_path: f64, wavelength: f64) -> f64 {
    if delta_path <= 0.0 {
        return 0.0;
    }
    let n = delta_path / wavelength;
    if n < 1e-10 {
        return 0.0;
    }
    -10.0 * (2.0 * PI * n).log10()
}
/// Powell-Howe vortex sound source term (2D): L = ρ (ω × u) · ê_r / (2π c₀²r)
///
/// `omega_z` is vorticity (z-component), `u` is velocity, `r_hat` is observer
/// direction, `r` is distance.
pub fn vortex_sound_pressure(
    rho: f64,
    omega_z: f64,
    u: [f64; 2],
    r_hat: [f64; 2],
    r: f64,
    c0: f64,
) -> f64 {
    let cross_x = -omega_z * u[1];
    let cross_y = omega_z * u[0];
    let dot = cross_x * r_hat[0] + cross_y * r_hat[1];
    rho * dot / (2.0 * PI * c0 * c0 * r)
}
/// Howe's acoustic analogy: mean-square acoustic pressure from turbulent vorticity.
///
/// p²_rms ≈ ρ² * U⁴ * Ω² * L³ / (c₀² * r²)
/// where U is velocity scale, Ω is vorticity scale, L is length scale.
pub fn howe_vortex_sound_power(rho: f64, u: f64, omega: f64, length: f64, c0: f64) -> f64 {
    let u2 = u * u;
    let u4 = u2 * u2;
    rho * rho * u4 * omega * omega * length * length * length / (c0 * c0)
}
/// Acoustic pressure from a point vortex convecting past an observer.
///
/// p ≈ ρ Γ Ω / (2π c₀ r) * cos(θ - θ_v)
/// where Γ is circulation, Ω is angular frequency, θ is observer angle.
pub fn point_vortex_acoustic_pressure(
    rho: f64,
    circulation: f64,
    omega: f64,
    r: f64,
    theta: f64,
    theta_v: f64,
    c0: f64,
) -> f64 {
    rho * circulation * omega / (2.0 * PI * c0 * r) * (theta - theta_v).cos()
}
/// Kármán vortex street shedding frequency (same as Strouhal relation).
pub fn karman_shedding_frequency(u: f64, d: f64, st: f64) -> f64 {
    st * u / d
}
/// Radiation resistance for a baffled piston of radius `a` in the low-frequency limit.
///
/// R_rad ≈ ρ c₀ π a² (k a)²  (k a << 1)
pub fn piston_radiation_resistance_low_freq(rho: f64, c0: f64, a: f64, k: f64) -> f64 {
    let ka = k * a;
    rho * c0 * PI * a * a * ka * ka
}
/// Radiation reactance for a baffled piston in the low-frequency limit.
///
/// X_rad ≈ ρ c₀ π a² (8/(3π)) * k a
pub fn piston_radiation_reactance_low_freq(rho: f64, c0: f64, a: f64, k: f64) -> f64 {
    let ka = k * a;
    rho * c0 * PI * a * a * (8.0 / (3.0 * PI)) * ka
}
/// Radiation efficiency of a planar source: σ = W / (ρ c₀ S `v²`).
///
/// Returns σ given radiated power `w`, surface area `s`, velocity rms `v_rms`.
pub fn planar_radiation_efficiency(w: f64, rho: f64, c0: f64, s: f64, v_rms: f64) -> f64 {
    let denominator = rho * c0 * s * v_rms * v_rms;
    if denominator.abs() < 1e-30 {
        return 0.0;
    }
    w / denominator
}
/// Modal overlap factor M = n(f) * Δf_loss / Δf_band.
///
/// `modal_density` is modes per Hz, `loss_factor` is η (dimensionless damping),
/// `freq` is centre frequency.
pub fn modal_overlap_factor(modal_density: f64, loss_factor: f64, freq: f64) -> f64 {
    PI * freq * loss_factor * modal_density
}
/// Room constant R = S α / (1 - α) where S is total surface area, α is mean absorption.
pub fn room_constant(s: f64, alpha: f64) -> f64 {
    if (1.0 - alpha).abs() < 1e-14 {
        return f64::INFINITY;
    }
    s * alpha / (1.0 - alpha)
}
/// Extract fluctuating density from LBM distribution functions.
///
/// ρ' = Σ_i f_i - ρ₀
pub fn lbm_density_fluctuation(f: &[f64; 9], rho0: f64) -> f64 {
    f.iter().sum::<f64>() - rho0
}
/// Compute acoustic pressure fluctuation from LBM: p' = cs² * ρ'
pub fn lbm_acoustic_pressure(f: &[f64; 9], rho0: f64, cs2: f64) -> f64 {
    cs2 * lbm_density_fluctuation(f, rho0)
}
/// Compute non-equilibrium stress tensor component P_xy from LBM distributions.
///
/// P_xy^{neq} = Σ_i (f_i - f_i^eq) * c_ix * c_iy
pub fn lbm_neq_stress_xy(f: &[f64; 9], f_eq: &[f64; 9], velocities: &[[i32; 2]; 9]) -> f64 {
    f.iter()
        .zip(f_eq.iter())
        .zip(velocities.iter())
        .map(|((fi, feqi), ci)| (fi - feqi) * ci[0] as f64 * ci[1] as f64)
        .sum()
}
/// Compute the acoustic Courant–Friedrichs–Lewy number for LBM.
///
/// CFL_ac = c₀ * Δt / Δx
pub fn acoustic_cfl(c0: f64, dt: f64, dx: f64) -> f64 {
    c0 * dt / dx
}
/// Acoustic absorption length l_abs = 1 / (2 * α_att)
/// where α_att is the spatial attenuation coefficient (Np/m).
pub fn acoustic_absorption_length(alpha_att: f64) -> f64 {
    if alpha_att.abs() < 1e-30 {
        return f64::INFINITY;
    }
    1.0 / (2.0 * alpha_att)
}
/// LBM acoustic relaxation parameter: τ_ac = 0.5 + c0 / (cs * c_lat)
///
/// Approximate formula relating physical sound speed `c0`, lattice sound speed
/// `cs`, and characteristic lattice velocity `c_lat`.
pub fn lbm_acoustic_tau(c0: f64, cs: f64, c_lat: f64) -> f64 {
    0.5 + c0 / (cs * c_lat)
}
/// Geometric spreading loss for a spherical wave: ΔSpl = -20*log10(r2/r1).
pub fn spherical_spreading_loss(r1: f64, r2: f64) -> f64 {
    if r1 <= 0.0 || r2 <= 0.0 {
        return 0.0;
    }
    -20.0 * (r2 / r1).log10()
}
/// Cylindrical spreading loss: ΔSpl = -10*log10(r2/r1).
pub fn cylindrical_spreading_loss(r1: f64, r2: f64) -> f64 {
    if r1 <= 0.0 || r2 <= 0.0 {
        return 0.0;
    }
    -10.0 * (r2 / r1).log10()
}
/// Excess attenuation due to ground reflection (simplified coherent model).
///
/// Uses path-length difference: Δ = r_reflected - r_direct
pub fn ground_reflection_excess_attenuation(
    p_direct: f64,
    path_diff: f64,
    freq: f64,
    c0: f64,
    ground_impedance: f64,
) -> f64 {
    let k = 2.0 * PI * freq / c0;
    let r_coeff = (ground_impedance - 1.0) / (ground_impedance + 1.0);
    let p_reflected = r_coeff * p_direct * (-(k * path_diff).powi(2)).exp();
    let p_total = p_direct + p_reflected;
    20.0 * (p_total.abs() / p_direct.abs()).log10()
}
/// Sutherland–Bass atmospheric absorption coefficient (m⁻¹) — simplified.
///
/// Uses the ISO 9613-1 approximation for air at temperature T_K (Kelvin) and
/// relative humidity rh (0–100).  Returns the energy attenuation in dB/m.
pub fn atmospheric_absorption_db_per_m(freq: f64, t_k: f64, rh: f64) -> f64 {
    let t_rel = t_k / 293.15;
    let h = rh * (t_rel).powf(-3.5) * (-2239.1 / t_k).exp().max(1e-30);
    let f_ro = 24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h);
    let f_rn =
        t_rel.powf(-0.5) * (9.0 + 280.0 * h * (-4.17 * (t_rel.powf(-1.0 / 3.0) - 1.0)).exp());
    let alpha = 1.84e-11 / t_rel.sqrt() * (freq / 1000.0).powi(2)
        + t_rel.powf(-2.5)
            * (0.01275 * (-2239.1 / t_k).exp() * (f_ro + freq * freq / f_ro).recip()
                + 0.1068 * (-3352.0 / t_k).exp() * (f_rn + freq * freq / f_rn).recip())
            * (freq / 1000.0).powi(2);
    alpha * 8.686
}
/// Add incoherent SPL contributions from multiple sources (energy summation).
pub fn add_incoherent_spl(spls: &[f64]) -> f64 {
    if spls.is_empty() {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = spls.iter().map(|&s| 10.0_f64.powf(s / 10.0)).sum();
    10.0 * sum.log10()
}
#[cfg(test)]
mod tests_fwh_and_noise {
    use super::*;

    #[test]
    fn test_fwh_thickness_noise_normal_velocity() {
        let panel = FwhSurfacePanel {
            area: 1.0,
            normal: [1.0, 0.0, 0.0],
            velocity: [10.0, 0.0, 0.0],
            delta_p: 0.0,
            r_hat: [1.0, 0.0, 0.0],
            r: 1.0,
        };
        let p = fwh_thickness_noise(&panel, 1.2, 340.0);
        assert!(p.abs() > 0.0, "FWH thickness noise should be non-zero: {p}");
    }
    #[test]
    fn test_fwh_thickness_noise_zero_velocity() {
        let panel = FwhSurfacePanel {
            area: 1.0,
            normal: [1.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
            delta_p: 100.0,
            r_hat: [1.0, 0.0, 0.0],
            r: 1.0,
        };
        let p = fwh_thickness_noise(&panel, 1.2, 340.0);
        assert_eq!(p, 0.0, "Zero velocity → zero thickness noise");
    }
    #[test]
    fn test_fwh_loading_noise_nonzero() {
        let panel = FwhSurfacePanel {
            area: 1.0,
            normal: [1.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
            delta_p: 100.0,
            r_hat: [1.0, 0.0, 0.0],
            r: 1.0,
        };
        let p = fwh_loading_noise(&panel, 340.0);
        assert!(p.abs() > 0.0, "FWH loading noise: {p}");
    }
    #[test]
    fn test_fwh_total_pressure_additive() {
        let panel = FwhSurfacePanel {
            area: 1.0,
            normal: [1.0, 0.0, 0.0],
            velocity: [1.0, 0.0, 0.0],
            delta_p: 1.0,
            r_hat: [1.0, 0.0, 0.0],
            r: 1.0,
        };
        let t = fwh_thickness_noise(&panel, 1.2, 340.0);
        let l = fwh_loading_noise(&panel, 340.0);
        let panel2 = FwhSurfacePanel {
            area: 1.0,
            normal: [1.0, 0.0, 0.0],
            velocity: [1.0, 0.0, 0.0],
            delta_p: 1.0,
            r_hat: [1.0, 0.0, 0.0],
            r: 1.0,
        };
        let total = fwh_total_pressure(&[panel2], 1.2, 340.0);
        assert!(
            (total - (t + l)).abs() < 1e-12,
            "Total pressure mismatch: {total} vs {}",
            t + l
        );
    }
    #[test]
    fn test_near_field_monopole_zero_distance() {
        let p = near_field_monopole_pressure(1.0, 1.0, 1.0, 0.0, 340.0);
        assert_eq!(p, 0.0);
    }
    #[test]
    fn test_near_field_monopole_large_r_dominates() {
        let p1 = near_field_monopole_pressure(1.0, 100.0, 0.0, 1.0, 340.0);
        let p2 = near_field_monopole_pressure(1.0, 100.0, 0.0, 10.0, 340.0);
        assert!(p1.abs() > p2.abs(), "Near-field: p1={p1}, p2={p2}");
    }
    #[test]
    fn test_near_field_dipole_cos_zero() {
        let p = near_field_dipole_pressure(1.0, 1.0, 1.0, 0.0, 340.0);
        assert_eq!(p, 0.0, "cos_theta=0 → zero dipole pressure");
    }
    #[test]
    fn test_trailing_edge_noise_peak_at_st_02() {
        let spl_02 = trailing_edge_noise_spl(0.2 * 10.0 / 0.01, 10.0, 0.01, 1.0, 1.2, 340.0);
        let spl_10 = trailing_edge_noise_spl(1.0 * 10.0 / 0.01, 10.0, 0.01, 1.0, 1.2, 340.0);
        assert!(
            spl_02 > spl_10,
            "Peak near St=0.2: spl_02={spl_02}, spl_10={spl_10}"
        );
    }
    #[test]
    fn test_trailing_edge_noise_zero_velocity() {
        let spl = trailing_edge_noise_spl(100.0, 0.0, 0.01, 1.0, 1.2, 340.0);
        assert!(spl.is_infinite() || spl < 0.0);
    }
    #[test]
    fn test_a_weighting_1khz_near_zero() {
        let aw = a_weighting_db(1000.0);
        assert!(
            (aw - 0.0).abs() < 1.0,
            "A-weighting at 1kHz should be ~0 dB: {aw}"
        );
    }
    #[test]
    fn test_a_weighting_low_freq_attenuated() {
        let aw_low = a_weighting_db(50.0);
        let aw_mid = a_weighting_db(1000.0);
        assert!(
            aw_low < aw_mid,
            "Low freq A-weighting should be lower: {aw_low} vs {aw_mid}"
        );
    }
    #[test]
    fn test_overall_a_weighted_spl_empty() {
        let val = overall_a_weighted_spl(&[], &[]);
        assert!(val.is_infinite());
    }
    #[test]
    fn test_overall_a_weighted_spl_positive() {
        let freqs = [1000.0, 2000.0, 4000.0];
        let spls = [60.0, 58.0, 55.0];
        let oaspl = overall_a_weighted_spl(&freqs, &spls);
        assert!(oaspl > 50.0, "OASPL should be above 50 dB: {oaspl}");
    }
    #[test]
    fn test_aeolian_tone_freq_strouhal_strouhal() {
        let f = aeolian_tone_freq_strouhal(10.0, 0.05, 0.2);
        assert!((f - 40.0).abs() < 1e-10, "f={f}");
    }
    #[test]
    fn test_cylinder_strouhal_roundtrip() {
        let st = cylinder_strouhal(40.0, 10.0, 0.05);
        assert!((st - 0.2).abs() < 1e-10, "st={st}");
    }
    #[test]
    fn test_tonal_lift_noise_spl_finite() {
        let spl = tonal_lift_noise_spl(100.0, 2.0 * PI * 200.0, 1.0, 1.2, 340.0);
        assert!(spl.is_finite(), "Tonal noise SPL should be finite: {spl}");
    }
    #[test]
    fn test_lock_in_bandwidth_width() {
        let (lo, hi) = lock_in_bandwidth(0.2);
        assert!(hi > lo, "lock-in: lo={lo}, hi={hi}");
        let center = 1.0 / 0.2;
        assert!((lo + hi) / 2.0 - center < 0.01, "Lock-in center mismatch");
    }
    #[test]
    fn test_acoustic_intensity_rms_plane_wave() {
        let p_rms = 1.0;
        let rho0 = 1.2;
        let c0 = 340.0;
        let i = acoustic_intensity_rms(p_rms, rho0, c0);
        let expected = 1.0 / (2.0 * rho0 * c0);
        assert!((i - expected).abs() < 1e-12, "I={i} expected={expected}");
    }
    #[test]
    fn test_acoustic_intensity_rms_vector_length() {
        let iv = acoustic_intensity_rms_vector(1.0, 1.2, 340.0, [1.0, 0.0, 0.0]);
        let mag = (iv[0] * iv[0] + iv[1] * iv[1] + iv[2] * iv[2]).sqrt();
        let i_scalar = acoustic_intensity_rms(1.0, 1.2, 340.0);
        assert!((mag - i_scalar).abs() < 1e-12, "mag={mag}");
    }
    #[test]
    fn test_acoustic_energy_density_rms() {
        let e = acoustic_energy_density_rms(1.0, 1.2, 340.0);
        assert!(e > 0.0, "E={e}");
    }
    #[test]
    fn test_acoustic_power_unit_area() {
        let i = acoustic_intensity_rms(1.0, 1.2, 340.0);
        let w = acoustic_power(1.0, 1.2, 340.0, 1.0);
        assert!((w - i).abs() < 1e-12, "W={w}, I={i}");
    }
    #[test]
    fn test_sound_power_level_db_reference() {
        let lw = sound_power_level_db(1e-12);
        assert!(
            (lw - 0.0).abs() < 1e-10,
            "LW at ref power should be 0 dB: {lw}"
        );
    }
    #[test]
    fn test_sound_power_level_db_zero_power() {
        let lw = sound_power_level_db(0.0);
        assert!(lw.is_infinite());
    }
    #[test]
    fn test_sound_power_level_db_roundtrip() {
        let w = 1e-6;
        let lw = sound_power_level_db(w);
        let w2 = sound_power_from_level(lw);
        assert!((w2 - w).abs() < 1e-18, "Roundtrip: w={w}, w2={w2}");
    }
    #[test]
    fn test_directivity_index_omnidirectional() {
        let total_power = 1.0;
        let r = 1.0;
        let i_omni = total_power / (4.0 * PI * r * r);
        let di = directivity_index(i_omni, total_power, r);
        assert!(di.abs() < 1e-10, "Omnidirectional DI should be 0: {di}");
    }
    #[test]
    fn test_sphere_radiated_power_positive() {
        let w = sphere_radiated_power(0.01, 2.0 * PI * 1000.0, 0.05, 1.2, 340.0);
        assert!(w > 0.0, "Sphere radiated power: {w}");
    }
    #[test]
    fn test_sphere_radiated_power_scales_with_velocity() {
        let w1 = sphere_radiated_power(0.01, 2.0 * PI * 1000.0, 0.05, 1.2, 340.0);
        let w2 = sphere_radiated_power(0.02, 2.0 * PI * 1000.0, 0.05, 1.2, 340.0);
        assert!((w2 / w1 - 4.0).abs() < 1e-10, "Power ratio: {}", w2 / w1);
    }
}
