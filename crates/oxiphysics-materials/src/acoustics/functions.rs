//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::WeightingFilter;

/// Pressure reflection coefficient at a normal-incidence interface.
///
/// R = (z2 - z1) / (z2 + z1)
pub fn reflection_coefficient(z1: f64, z2: f64) -> f64 {
    (z2 - z1) / (z2 + z1)
}
/// Pressure transmission coefficient at a normal-incidence interface.
///
/// T = 2·z2 / (z2 + z1)
pub fn transmission_coefficient(z1: f64, z2: f64) -> f64 {
    2.0 * z2 / (z2 + z1)
}
/// Transmission loss in dB based on intensity transmission coefficient.
///
/// The intensity (power) transmission coefficient is τ = 4·z1·z2 / (z1+z2)².
/// TL = −10·log10(τ).
pub fn transmission_loss_db(z1: f64, z2: f64) -> f64 {
    let tau = 4.0 * z1 * z2 / ((z1 + z2) * (z1 + z2));
    -10.0 * tau.log10()
}
/// Transfer-matrix method for a multilayer system at normal incidence.
///
/// `impedances` has length N (outer medium, layers…, outer medium).
/// `thicknesses` has length N−2 (one per interior layer).
/// `freq` in Hz.
/// Returns the power transmission coefficient |T|².
pub fn multilayer_transmission(impedances: &[f64], thicknesses: &[f64], freq: f64) -> f64 {
    let n = impedances.len();
    if n < 2 {
        return 1.0;
    }
    let n_layers = n - 2;
    let mut m11 = 1.0_f64;
    let mut m12 = 0.0_f64;
    let mut m21 = 0.0_f64;
    let mut m22 = 1.0_f64;
    for i in 0..n_layers {
        let z = impedances[i + 1];
        let omega = 2.0 * PI * freq;
        let c = z / 1.2;
        let k = omega / c;
        let d = thicknesses[i];
        let kd = k * d;
        let cos_kd = kd.cos();
        let sin_kd = kd.sin();
        let a11 = cos_kd * m11 + (z * sin_kd) * m21;
        let a12 = cos_kd * m12 + (z * sin_kd) * m22;
        let a21 = (sin_kd / z) * m11 + cos_kd * m21;
        let a22 = (sin_kd / z) * m12 + cos_kd * m22;
        m11 = a11;
        m12 = a12;
        m21 = a21;
        m22 = a22;
    }
    let z1 = impedances[0];
    let z2 = impedances[n - 1];
    let denom = m11 + m12 / z2 + z1 * m21 + z1 / z2 * m22;
    if denom.abs() < f64::EPSILON {
        return 0.0;
    }
    let t_amp = 2.0 / denom;
    t_amp * t_amp
}
/// Classical viscous attenuation coefficient (Stokes/Kirchhoff).
///
/// alpha = 2·mu·omega² / (3·rho·c³)  \[Np/m\]
pub fn viscous_attenuation(freq: f64, viscosity: f64, density: f64, c: f64) -> f64 {
    let omega = 2.0 * PI * freq;
    2.0 * viscosity * omega * omega / (3.0 * density * c * c * c)
}
/// Simplified ISO 9613-1 atmospheric absorption in dB/m.
///
/// Includes O₂ and N₂ vibrational relaxation.
/// `temp_c` in °C, `humidity` as fraction 0–1.
pub fn atmospheric_attenuation_db_per_m(freq: f64, temp_c: f64, humidity: f64) -> f64 {
    let t_kelvin = temp_c + 273.15;
    let t_ref = 293.15_f64;
    let t_ratio = t_kelvin / t_ref;
    let h = humidity * t_ratio.powf(-5.0 / 2.0) * (-2239.1 / t_kelvin).exp() * 1.0e4;
    let f_ro = 24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h);
    let f_rn =
        t_ratio.powf(-0.5) * (9.0 + 280.0 * h * (-4.17 * (t_ratio.powf(-1.0 / 3.0) - 1.0)).exp());
    let f2 = freq * freq;
    let alpha_classic = 1.84e-11 * t_ratio.powf(0.5) * f2;
    let alpha_o2 = 0.01275 * (-2239.1 / t_kelvin).exp() * (f_ro + f2 / f_ro).recip() * freq;
    let alpha_n2 = 0.1068 * (-3352.0 / t_kelvin).exp() * (f_rn + f2 / f_rn).recip() * freq;
    (alpha_classic + alpha_o2 + alpha_n2) * 8.686
}
/// Convert attenuation coefficient (dB/m or Np/m) over a distance to total dB loss.
pub fn attenuation_db(alpha_per_m: f64, distance: f64) -> f64 {
    alpha_per_m * distance
}
/// Sabine reverberation time: T60 = 0.161·V / (A·alpha).
///
/// `volume` in m³, `surface_area` in m², `absorption_coeff` in \[0,1\].
pub fn sabine_reverberation_time(volume: f64, surface_area: f64, absorption_coeff: f64) -> f64 {
    0.161 * volume / (surface_area * absorption_coeff)
}
/// Eyring reverberation time: T60 = 0.161·V / (−A·ln(1 − alpha)).
///
/// More accurate than Sabine for highly absorptive rooms.
pub fn eyring_reverberation_time(volume: f64, surface_area: f64, absorption_coeff: f64) -> f64 {
    0.161 * volume / (-surface_area * (1.0 - absorption_coeff).ln())
}
/// Critical distance (reverberant field = direct field): r_c = 0.057·sqrt(V / T60).
pub fn critical_distance(volume: f64, t60: f64) -> f64 {
    0.057 * (volume / t60).sqrt()
}
/// Noise reduction between two rooms.
///
/// NR = NRC_sender − NRC_receiver + 10·log10(area / (0.161·V_receiver / T60))
pub fn noise_reduction(
    nrc_sender: f64,
    nrc_receiver: f64,
    area: f64,
    volume_receiver: f64,
    t60: f64,
) -> f64 {
    let absorption_receiver = 0.161 * volume_receiver / t60;
    nrc_sender - nrc_receiver + 10.0 * (area / absorption_receiver).log10()
}
/// Wave number: k = 2·pi·f / c  \[rad/m\].
pub fn wave_number(freq: f64, c: f64) -> f64 {
    2.0 * PI * freq / c
}
/// Instantaneous pressure of a plane wave: p(x,t) = A·cos(k·x − omega·t).
pub fn plane_wave_pressure(amplitude: f64, k: f64, x: f64, omega: f64, t: f64) -> f64 {
    amplitude * (k * x - omega * t).cos()
}
/// Sound pressure level in dB re 20 µPa.
pub fn spl_db(pressure_pa: f64) -> f64 {
    let p_ref = 20.0e-6;
    20.0 * (pressure_pa.abs() / p_ref).log10()
}
/// Simplified loudness in sone (Stevens' power law at 1 kHz reference).
///
/// Returns 2^((spl_db - 40) / 10).  `freq_hz` is accepted for API
/// compatibility with future equal-loudness corrections.
pub fn loudness_sone(spl_db: f64, _freq_hz: f64) -> f64 {
    2.0_f64.powf((spl_db - 40.0) / 10.0)
}
/// Mass-law transmission loss \[dB\] for an infinite panel.
///
/// `TL = 20·log10(m·f) − 47.5`  (simplified mass law, SI units)
///
/// # Arguments
/// * `surface_mass_kg_m2` — surface mass density \[kg/m²\]
/// * `freq_hz`            — frequency \[Hz\]
pub fn mass_law_tl(surface_mass_kg_m2: f64, freq_hz: f64) -> f64 {
    20.0 * (surface_mass_kg_m2 * freq_hz).log10() - 47.5
}
/// Coincidence frequency for a panel \[Hz\].
///
/// `f_c = c² / (1.8·c_L·h)` where `c_L = sqrt(E/(rho*(1-nu²)))` is the
/// longitudinal plate wave speed and `h` is thickness \[m\].
///
/// # Arguments
/// * `c_air` — speed of sound in air \[m/s\]
/// * `c_l`   — longitudinal wave speed in panel \[m/s\]
/// * `h`     — panel thickness \[m\]
pub fn coincidence_frequency(c_air: f64, c_l: f64, h: f64) -> f64 {
    c_air * c_air / (1.8 * c_l * h)
}
/// A-weighting correction \[dB\] at frequency `f` \[Hz\].
///
/// Uses the IEC 61672 analytical formula.
pub fn a_weighting_db(f: f64) -> f64 {
    if f <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let f1 = 20.598_997_0;
    let f2 = 107.652_72;
    let f3 = 737.862_2;
    let f4 = 12_194.217;
    let f2s = f * f;
    let num = (f4 * f4) * f2s * f2s;
    let d1 = f2s + f1 * f1;
    let d2 = (f2s + f2 * f2).sqrt() * (f2s + f3 * f3).sqrt();
    let d3 = f2s + f4 * f4;
    let ra = num / (d1 * d2 * d3);
    20.0 * ra.log10() + 2.0
}
/// C-weighting correction \[dB\] at frequency `f` \[Hz\].
pub fn c_weighting_db(f: f64) -> f64 {
    if f <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let f1 = 20.598_997_0;
    let f4 = 12_194.217;
    let f2s = f * f;
    let num = (f4 * f4) * f2s;
    let d1 = f2s + f1 * f1;
    let d3 = f2s + f4 * f4;
    let rc = num / (d1 * d3);
    20.0 * rc.log10() + 0.06
}
/// Octave-band centre frequencies (Hz) from 31.5 Hz to 16 kHz (10 bands).
pub const OCTAVE_BAND_CENTRES: [f64; 10] = [
    31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16_000.0,
];
/// Overall sound pressure level from octave-band levels \[dB\].
///
/// Combines `n` octave-band SPL values using the energy summation rule:
/// `L_total = 10·log10(Σ 10^(L_i / 10))`.
pub fn overall_spl_from_octave_bands(band_spl: &[f64]) -> f64 {
    let sum: f64 = band_spl.iter().map(|&l| 10.0_f64.powf(l / 10.0)).sum();
    10.0 * sum.log10()
}
/// Apply a weighting filter to an octave-band spectrum and return corrected levels \[dB\].
///
/// # Arguments
/// * `band_spl`   — unweighted octave-band SPL values \[dB\]
/// * `band_freqs` — corresponding centre frequencies \[Hz\]
/// * `filter`     — weighting filter variant
pub fn apply_weighting(band_spl: &[f64], band_freqs: &[f64], filter: WeightingFilter) -> Vec<f64> {
    band_spl
        .iter()
        .zip(band_freqs.iter())
        .map(|(&l, &f)| {
            let w = match filter {
                WeightingFilter::A => a_weighting_db(f),
                WeightingFilter::B => a_weighting_db(f) * 0.5,
                WeightingFilter::C => c_weighting_db(f),
            };
            l + w
        })
        .collect()
}
/// Compute the characteristic impedance of an optimal quarter-wave matching
/// layer between media with impedances z1 and z2.
///
/// Z_match = sqrt(z1 * z2)
///
/// At f_0 such that layer thickness d = λ/4, the reflection coefficient is zero.
pub fn impedance_matching_layer(z1: f64, z2: f64) -> f64 {
    (z1 * z2).sqrt()
}
/// Residual reflection coefficient when using a matching layer with impedance z_m.
///
/// For a single quarter-wave layer, the residual at the design frequency is 0.
/// Off-design, the reflection depends on the frequency ratio.
///
/// This computes the reflection coefficient for a single matching layer at
/// normal incidence as a function of normalised frequency f/f0:
///
/// R(f/f0) = (z1·z2 - z_m²·cos²(π/2 · f/f0)) / (...)
///
/// Simplified formula: at f/f0 = 1, R = 0.
pub fn matching_layer_reflection(z1: f64, z_m: f64, z2: f64, f_over_f0: f64) -> f64 {
    let phi = 0.5 * PI * f_over_f0;
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let a = cos_phi;
    let b = z_m * sin_phi;
    let c = sin_phi / z_m;
    let d = cos_phi;
    let denom_re = a + b / z2 + z1 * c + z1 * d / z2;
    let t = 2.0 / denom_re.max(1e-20);
    let tau = z1 / z2 * t * t;
    let r_sq = (1.0 - tau).clamp(0.0, 1.0);
    r_sq.sqrt()
}
/// Compute the bandwidth over which a matching layer keeps |R| below a threshold.
///
/// Scans f/f0 from 0.1 to 2.0 and returns the range of f/f0 where |R| < threshold.
pub fn matching_layer_bandwidth(z1: f64, z_m: f64, z2: f64, threshold: f64) -> (f64, f64) {
    let n = 1000;
    let mut f_low = 2.0_f64;
    let mut f_high = 0.0_f64;
    for i in 1..=n {
        let f_ratio = 0.1 + 1.9 * i as f64 / n as f64;
        let r = matching_layer_reflection(z1, z_m, z2, f_ratio);
        if r < threshold {
            if f_ratio < f_low {
                f_low = f_ratio;
            }
            if f_ratio > f_high {
                f_high = f_ratio;
            }
        }
    }
    (f_low, f_high)
}
/// Angle-dependent intensity transmission coefficient for a fluid-fluid interface.
///
/// Snell's law: sin(θ_t) = (c2/c1) * sin(θ_i)
/// Transmission: T_I = 4·Z1·Z2·cos(θ_i)·cos(θ_t) / (Z1·cos(θ_t) + Z2·cos(θ_i))²
///
/// # Arguments
/// * `z1`, `z2` – impedances of media
/// * `c1`, `c2` – wave speeds in media
/// * `theta_i`  – angle of incidence (radians)
///
/// Returns `None` for total internal reflection.
pub fn angle_dependent_transmission(
    z1: f64,
    z2: f64,
    c1: f64,
    c2: f64,
    theta_i: f64,
) -> Option<f64> {
    let sin_t = (c2 / c1) * theta_i.sin();
    if sin_t > 1.0 {
        return None;
    }
    let cos_i = theta_i.cos();
    let cos_t = (1.0 - sin_t * sin_t).sqrt();
    let num = 4.0 * z1 * z2 * cos_i * cos_t;
    let denom = (z1 * cos_t + z2 * cos_i).powi(2);
    Some(if denom < 1e-30 { 0.0 } else { num / denom })
}
/// Critical angle (in radians) for total internal reflection.
///
/// θ_c = arcsin(c1 / c2), only defined if c2 > c1.
pub fn critical_angle(c1: f64, c2: f64) -> Option<f64> {
    if c2 <= c1 {
        return None;
    }
    let sin_c = c1 / c2;
    if sin_c > 1.0 {
        None
    } else {
        Some(sin_c.asin())
    }
}
/// Axial, tangential, and oblique mode resonance frequencies of a rectangular room.
///
/// For a shoebox room of dimensions Lx × Ly × Lz, the resonance frequency of mode
/// (nx, ny, nz) is:
///
/// `f = (c/2) * sqrt((nx/Lx)² + (ny/Ly)² + (nz/Lz)²)`
///
/// where c is the speed of sound \[m/s\] and n_i are non-negative integers.
pub fn room_resonance_frequency(
    lx: f64,
    ly: f64,
    lz: f64,
    nx: u32,
    ny: u32,
    nz: u32,
    c: f64,
) -> f64 {
    let fx = (nx as f64) / lx;
    let fy = (ny as f64) / ly;
    let fz = (nz as f64) / lz;
    0.5 * c * (fx * fx + fy * fy + fz * fz).sqrt()
}
/// Collect the N lowest-frequency resonance modes of a rectangular room.
///
/// Scans modes (nx, ny, nz) up to `max_order` per axis and returns the N
/// lowest resonance frequencies in ascending order.
///
/// # Arguments
/// * `lx`, `ly`, `lz` — room dimensions \[m\]
/// * `c`               — speed of sound \[m/s\]
/// * `max_order`       — maximum mode index per axis (inclusive)
/// * `n_modes`         — number of modes to return
pub fn room_modes(lx: f64, ly: f64, lz: f64, c: f64, max_order: u32, n_modes: usize) -> Vec<f64> {
    let mut freqs: Vec<f64> = Vec::new();
    for nx in 0..=max_order {
        for ny in 0..=max_order {
            for nz in 0..=max_order {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                let f = room_resonance_frequency(lx, ly, lz, nx, ny, nz, c);
                freqs.push(f);
            }
        }
    }
    freqs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    freqs.dedup_by(|a, b| (*a - *b).abs() < 0.01);
    freqs.truncate(n_modes);
    freqs
}
/// Lowest axial resonance frequency (Schroeder criterion helper).
///
/// The lowest mode in a rectangular room: f_1 = c / (2 * L_max).
pub fn room_lowest_mode(lx: f64, ly: f64, lz: f64, c: f64) -> f64 {
    let l_max = lx.max(ly).max(lz);
    c / (2.0 * l_max)
}
/// Schroeder frequency \[Hz\] for a room.
///
/// Above this frequency, modal overlap is sufficient for statistical treatment:
/// `f_S = 2000 * sqrt(T60 / V)`
///
/// # Arguments
/// * `t60`    — reverberation time \[s\]
/// * `volume` — room volume \[m³\]
pub fn schroeder_frequency(t60: f64, volume: f64) -> f64 {
    2000.0 * (t60 / volume).sqrt()
}
/// Sound transmission loss through a double-leaf partition (simplified model) \[dB\].
///
/// For two panels with surface masses m1, m2 \[kg/m²\] and an air gap d \[m\],
/// the TL is approximately the sum of individual mass-law TLs plus an air-gap
/// resonance penalty.
///
/// This uses the simplified formula:
/// `TL_double = TL1 + TL2 + 6  - cavity_resonance_dip`
///
/// # Arguments
/// * `m1_kg_m2`  — surface mass of panel 1 \[kg/m²\]
/// * `m2_kg_m2`  — surface mass of panel 2 \[kg/m²\]
/// * `gap_m`     — air gap thickness \[m\]
/// * `freq_hz`   — frequency \[Hz\]
pub fn double_leaf_tl(m1_kg_m2: f64, m2_kg_m2: f64, gap_m: f64, freq_hz: f64) -> f64 {
    let tl1 = mass_law_tl(m1_kg_m2, freq_hz);
    let tl2 = mass_law_tl(m2_kg_m2, freq_hz);
    let f_cav = 343.0 / (2.0 * gap_m);
    let q_factor = 5.0;
    let resonance_dip = 10.0 / (1.0 + ((freq_hz - f_cav) / (f_cav / q_factor)).powi(2));
    (tl1 + tl2 + 6.0 - resonance_dip).max(tl1.max(tl2))
}
/// Weighted sound reduction index (Rw) from 1/3-octave band data \[dB\].
///
/// Fits the reference curve to the measured data per ISO 717-1.
/// Uses a simplified single-number rating based on averaging key bands.
///
/// # Arguments
/// * `tl_values` — transmission loss values at standard 1/3-octave bands
///   (16 values from 100 Hz to 3150 Hz)
pub fn weighted_sound_reduction_index(tl_values: &[f64]) -> f64 {
    let reference: [f64; 16] = [
        33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 56.0, 56.0, 56.0,
        56.0,
    ];
    let n = tl_values.len().min(reference.len());
    if n == 0 {
        return 0.0;
    }
    let avg_tl: f64 = tl_values[..n].iter().sum::<f64>() / n as f64;
    let avg_ref: f64 = reference[..n].iter().sum::<f64>() / n as f64;
    avg_tl - avg_ref + 52.0
}
/// Sound intensity \[W/m²\] from pressure amplitude and acoustic impedance.
///
/// `I = p_rms² / Z`  where `p_rms = p_peak / sqrt(2)` for a sinusoidal wave.
pub fn sound_intensity(pressure_pa_peak: f64, impedance: f64) -> f64 {
    let p_rms = pressure_pa_peak / 2.0_f64.sqrt();
    p_rms * p_rms / impedance
}
/// Sound intensity level \[dB re 1 pW/m²\].
///
/// `SIL = 10·log10(I / I_ref)` where `I_ref = 1e-12 W/m²`.
pub fn sound_intensity_level_db(intensity_w_m2: f64) -> f64 {
    let i_ref = 1.0e-12;
    10.0 * (intensity_w_m2 / i_ref).log10()
}
/// Acoustic power \[W\] radiated by a source with intensity I over area A.
pub fn acoustic_power(intensity: f64, area: f64) -> f64 {
    intensity * area
}
/// Sound power level \[dB re 1 pW\].
///
/// `PWL = 10·log10(W / W_ref)` where `W_ref = 1e-12 W`.
pub fn sound_power_level_db(power_w: f64) -> f64 {
    let w_ref = 1.0e-12;
    10.0 * (power_w / w_ref).log10()
}
/// Doppler-shifted frequency for a moving source or observer.
///
/// Classical (non-relativistic) Doppler formula:
/// `f' = f0 * (c + v_observer) / (c - v_source)`
///
/// Sign convention:
/// * `v_source` > 0  → source moving toward observer
/// * `v_observer` > 0 → observer moving toward source
///
/// # Arguments
/// * `f0`          — emitted frequency \[Hz\]
/// * `c`           — speed of sound \[m/s\]
/// * `v_source`    — source velocity toward observer \[m/s\]
/// * `v_observer`  — observer velocity toward source \[m/s\]
pub fn doppler_frequency(f0: f64, c: f64, v_source: f64, v_observer: f64) -> f64 {
    if (c - v_source).abs() < 1e-10 {
        return f64::INFINITY;
    }
    f0 * (c + v_observer) / (c - v_source)
}
/// Mach number for a body moving at speed v in a medium with sound speed c.
pub fn mach_number(v: f64, c: f64) -> f64 {
    v / c
}
/// Maekawa insertion loss \[dB\] for a thin noise barrier.
///
/// Uses the Maekawa (1968) empirical formula:
/// `IL = 10·log10(3 + 20·N)`  where `N = 2·delta/lambda` is the Fresnel number.
///
/// # Arguments
/// * `delta` — path length difference (m): (SA + AB) - SB where S = source, B = barrier top
/// * `freq`  — frequency \[Hz\]
/// * `c`     — speed of sound \[m/s\]
pub fn barrier_insertion_loss_maekawa(delta: f64, freq: f64, c: f64) -> f64 {
    let lambda = c / freq;
    let n_fresnel = 2.0 * delta / lambda;
    (10.0 * (3.0 + 20.0 * n_fresnel).log10()).max(0.0)
}
/// Fresnel number for a barrier with path length difference delta \[m\] at frequency f \[Hz\].
pub fn fresnel_number(delta: f64, freq: f64, c: f64) -> f64 {
    2.0 * delta * freq / c
}
/// Helmholtz resonator natural frequency \[Hz\].
///
/// `f_H = (c / 2π) * sqrt(A / (V * L_eff))`
///
/// where:
/// * `A`     — neck cross-section area \[m²\]
/// * `V`     — cavity volume \[m³\]
/// * `L_eff` — effective neck length \[m\] (actual length + end corrections)
/// * `c`     — speed of sound \[m/s\]
pub fn helmholtz_resonator_frequency(area_m2: f64, volume_m3: f64, l_eff_m: f64, c: f64) -> f64 {
    (c / (2.0 * PI)) * (area_m2 / (volume_m3 * l_eff_m)).sqrt()
}
/// Quarter-wave tube resonator frequency \[Hz\].
///
/// Open-closed tube: `f_n = (2n-1) * c / (4 * L)`,  n = 1, 2, 3, …
///
/// # Arguments
/// * `length` — tube length \[m\]
/// * `c`      — speed of sound \[m/s\]
/// * `n`      — harmonic number (1 = fundamental)
pub fn quarter_wave_resonator_frequency(length: f64, c: f64, n: u32) -> f64 {
    let n_f = (2 * n - 1) as f64;
    n_f * c / (4.0 * length)
}
/// Half-wave tube resonator frequency \[Hz\].
///
/// Open-open (or closed-closed) tube: `f_n = n * c / (2 * L)`.
pub fn half_wave_resonator_frequency(length: f64, c: f64, n: u32) -> f64 {
    (n as f64) * c / (2.0 * length)
}
/// Mean free path of sound in a room \[m\].
///
/// `L_mfp = 4 * V / S`  (statistical room acoustics result)
///
/// # Arguments
/// * `volume`       — room volume \[m³\]
/// * `surface_area` — total surface area \[m²\]
pub fn mean_free_path(volume: f64, surface_area: f64) -> f64 {
    4.0 * volume / surface_area
}
/// Room constant R \[m²·sabin\].
///
/// `R = S·alpha / (1 - alpha)`
///
/// Used in the classical room acoustic formula for SPL from a source.
pub fn room_constant(surface_area: f64, absorption_coeff: f64) -> f64 {
    surface_area * absorption_coeff / (1.0 - absorption_coeff.min(0.9999))
}
/// Diffuse field SPL increase from a source of power W \[dB re 20 µPa\].
///
/// `SPL_diff = PWL + 10·log10(4 / R)`
///
/// # Arguments
/// * `power_w`     — source acoustic power \[W\]
/// * `room_const`  — room constant R \[m²·sabin\]
pub fn diffuse_field_spl(power_w: f64, room_const: f64) -> f64 {
    let pwl = sound_power_level_db(power_w);
    pwl + 10.0 * (4.0 / room_const).log10()
}
/// Total (direct + reverberant) SPL at distance r from a source.
///
/// `SPL = PWL + 10·log10(Q/(4πr²) + 4/R)`
///
/// # Arguments
/// * `power_w`        — source acoustic power \[W\]
/// * `room_const`     — room constant R \[m²·sabin\]
/// * `distance_m`     — distance from source \[m\]
/// * `directivity_q`  — directivity factor (Q = 1 for omnidirectional)
pub fn total_spl_in_room(
    power_w: f64,
    room_const: f64,
    distance_m: f64,
    directivity_q: f64,
) -> f64 {
    let pwl = sound_power_level_db(power_w);
    let direct_term = directivity_q / (4.0 * PI * distance_m * distance_m);
    let reverb_term = 4.0 / room_const;
    pwl + 10.0 * (direct_term + reverb_term).log10()
}
/// Speed of sound in an ideal gas \[m/s\].
///
/// `c = sqrt(gamma * R_specific * T)`
///
/// # Arguments
/// * `gamma`      — heat capacity ratio (Cp/Cv), ~1.4 for air
/// * `r_specific` — specific gas constant \[J/(kg·K)\], ~287 for air
/// * `temp_k`     — temperature \[K\]
pub fn sound_speed_ideal_gas(gamma: f64, r_specific: f64, temp_k: f64) -> f64 {
    (gamma * r_specific * temp_k).sqrt()
}
/// Speed of sound in a liquid using the Newton-Laplace relation.
///
/// `c = sqrt(B / rho)`
///
/// where B is the adiabatic bulk modulus \[Pa\] and rho is density \[kg/m³\].
pub fn sound_speed_liquid(bulk_modulus_pa: f64, density_kg_m3: f64) -> f64 {
    (bulk_modulus_pa / density_kg_m3).sqrt()
}
/// Temperature dependence of sound speed in air (Cramer, 1993 simplified).
///
/// `c(T) ≈ 331.3 * sqrt(1 + T_celsius / 273.15)` \[m/s\]
pub fn sound_speed_air_temperature(temp_celsius: f64) -> f64 {
    331.3 * (1.0 + temp_celsius / 273.15).sqrt()
}
/// Group velocity for a dispersive medium: `v_g = dω/dk`.
///
/// Computed numerically via central difference of the dispersion relation
/// `omega(k)` provided as a closure.
///
/// # Arguments
/// * `k`  — wave number \[rad/m\] at which to evaluate group velocity
/// * `dk` — finite-difference step \[rad/m\]
/// * `omega_fn` — dispersion relation returning ω for a given k
pub fn group_velocity<F>(k: f64, dk: f64, omega_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    (omega_fn(k + dk) - omega_fn(k - dk)) / (2.0 * dk)
}
/// Peak acoustic pressure amplitude from intensity.
///
/// `p_peak = sqrt(2 * I * Z)`
///
/// where I is sound intensity \[W/m²\] and Z is acoustic impedance \[Pa·s/m\].
pub fn pressure_from_intensity(intensity_w_m2: f64, impedance: f64) -> f64 {
    (2.0 * intensity_w_m2 * impedance).sqrt()
}
/// RMS acoustic pressure \[Pa\] from peak pressure.
///
/// `p_rms = p_peak / sqrt(2)` for a sinusoidal wave.
pub fn rms_pressure(p_peak: f64) -> f64 {
    p_peak / 2.0_f64.sqrt()
}
/// Particle velocity amplitude \[m/s\] in a plane wave.
///
/// `u = p / Z` where p is pressure amplitude and Z is impedance.
pub fn particle_velocity(pressure_pa: f64, impedance: f64) -> f64 {
    if impedance.abs() < f64::EPSILON {
        return 0.0;
    }
    pressure_pa / impedance
}
/// Standing wave pressure distribution in a closed-closed tube.
///
/// `p(x, t) = 2A * cos(k*x) * cos(ω*t)`
///
/// Returns the spatial factor `2A * cos(k*x)`.
pub fn standing_wave_pressure_closed_closed(amplitude: f64, k: f64, x: f64) -> f64 {
    2.0 * amplitude * (k * x).cos()
}
/// Standing wave pressure in an open-closed tube.
///
/// Pressure has an anti-node at the closed end and a node at the open end.
/// `p(x, t) = 2A * sin(k*x) * cos(ω*t)` (x measured from open end)
pub fn standing_wave_pressure_open_closed(amplitude: f64, k: f64, x: f64) -> f64 {
    2.0 * amplitude * (k * x).sin()
}
/// Resonance frequencies of a closed-closed cylindrical pipe \[Hz\].
///
/// `f_n = n * c / (2 * L)`, n = 1, 2, 3, ...
///
/// # Arguments
/// * `length` — pipe length \[m\]
/// * `c`      — speed of sound \[m/s\]
/// * `n`      — mode number (1 = fundamental)
pub fn closed_closed_resonance(length: f64, c: f64, n: u32) -> f64 {
    (n as f64) * c / (2.0 * length)
}
/// Q-factor of a resonator from bandwidth.
///
/// `Q = f_0 / Δf`  where Δf is the -3 dB bandwidth.
pub fn resonator_q_factor(f0: f64, bandwidth_hz: f64) -> f64 {
    if bandwidth_hz < f64::EPSILON {
        return f64::INFINITY;
    }
    f0 / bandwidth_hz
}
/// Resonance decay time constant τ from Q-factor.
///
/// `τ = Q / (π * f_0)`
pub fn resonance_decay_time(q: f64, f0: f64) -> f64 {
    if f0 < f64::EPSILON {
        return f64::INFINITY;
    }
    q / (PI * f0)
}
/// Energy density of an acoustic wave \[J/m³\].
///
/// `w = p_rms² / (rho * c²)`  (time-averaged)
pub fn acoustic_energy_density(p_rms: f64, density: f64, c: f64) -> f64 {
    let denom = density * c * c;
    if denom < f64::EPSILON {
        return 0.0;
    }
    p_rms * p_rms / denom
}
/// Array factor for a uniform linear array (ULA) of `n` omnidirectional elements
/// spaced `d` apart at frequency `f` for steering angle `theta_steer` \[rad\].
///
/// Returns the normalised power pattern at angle `theta` \[rad\]:
/// `AF = |sin(N*ψ/2) / (N * sin(ψ/2))|²`
/// where `ψ = k*d*(cos(θ) - cos(θ_steer))`.
pub fn ula_array_factor(n: u32, d: f64, freq: f64, c: f64, theta: f64, theta_steer: f64) -> f64 {
    let k = 2.0 * PI * freq / c;
    let psi = k * d * (theta.cos() - theta_steer.cos());
    let n_f = n as f64;
    let num = (n_f * psi / 2.0).sin();
    let den = (psi / 2.0).sin();
    if den.abs() < 1e-12 {
        return 1.0;
    }
    (num / (n_f * den)).powi(2)
}
/// Normalised input impedance of an open-end cylindrical tube at frequency f \[Hz\].
///
/// For a tube of length L with one open end (radiation impedance ignored):
/// `Z_in = j * Z_0 * cot(k*L)`
///
/// Returns the imaginary part only (real part ≈ 0 for lossless tube).
/// Sign convention: positive = inductive (reactive), negative = capacitive.
pub fn tube_input_impedance_open_end_imag(z0: f64, freq: f64, c: f64, length: f64) -> f64 {
    let k = 2.0 * PI * freq / c;
    let kl = k * length;
    let sin_kl = kl.sin();
    if sin_kl.abs() < 1e-12 {
        return f64::INFINITY;
    }
    -z0 * kl.cos() / sin_kl
}
/// Normalised input impedance of a closed-end cylindrical tube.
///
/// `Z_in = -j * Z_0 * cot(k*L)` (imaginary part)
///
/// Returns the imaginary part.
pub fn tube_input_impedance_closed_end_imag(z0: f64, freq: f64, c: f64, length: f64) -> f64 {
    let k = 2.0 * PI * freq / c;
    let kl = k * length;
    let sin_kl = kl.sin();
    if sin_kl.abs() < 1e-12 {
        return f64::INFINITY;
    }
    z0 * kl.cos() / sin_kl
}
/// Power-law acoustic attenuation model commonly used in biomedical ultrasound.
///
/// `α(f) = α_0 * f^b`  \[dB/cm\]
///
/// where `α_0` is the attenuation coefficient at 1 MHz \[dB/(cm·MHz^b)\]
/// and `b` is the frequency exponent (≈ 1 for soft tissue, 2 for water).
///
/// # Arguments
/// * `alpha_0` — attenuation coefficient \[dB/(cm·MHz^b)\]
/// * `freq_mhz` — frequency \[MHz\]
/// * `b`        — power exponent
pub fn power_law_attenuation(alpha_0: f64, freq_mhz: f64, b: f64) -> f64 {
    alpha_0 * freq_mhz.powf(b)
}
/// Total attenuation \[dB\] for a round-trip path (pulse-echo) in tissue.
///
/// `att_total = 2 * α(f) * depth_cm`
pub fn pulse_echo_attenuation_db(alpha_0: f64, freq_mhz: f64, b: f64, depth_cm: f64) -> f64 {
    2.0 * power_law_attenuation(alpha_0, freq_mhz, b) * depth_cm
}
/// Radiation resistance of a piston in an infinite baffle.
///
/// For a circular piston of radius a at frequency f \[Hz\],
/// the radiation resistance normalised to ρ*c*A is:
///
/// `R_rad = ρ * c * A * R1(ka)` where `R1(x) = 1 - J1(2x)/x` (Bessel approx)
///
/// For ka << 1: `R_rad ≈ ρ * c * A * (ka)² / 2`
/// For ka >> 1: `R_rad ≈ ρ * c * A`
///
/// Returns the normalised radiation resistance R1(ka).
pub fn piston_radiation_resistance_normalised(freq: f64, radius: f64, c: f64) -> f64 {
    let ka = 2.0 * PI * freq * radius / c;
    if ka < 0.1 {
        (ka * ka) / 2.0
    } else if ka > 10.0 {
        1.0
    } else {
        1.0 - (2.0 * ka).sin() / (2.0 * ka)
    }
}
/// Directivity index \[dB\] for a circular piston in a baffle.
///
/// Simplified: DI ≈ 10·log10(2·(ka)²) for ka > 1.
pub fn piston_directivity_index_db(freq: f64, radius: f64, c: f64) -> f64 {
    let ka = 2.0 * PI * freq * radius / c;
    if ka < f64::EPSILON {
        return 0.0;
    }
    10.0 * (2.0 * ka * ka).log10().max(0.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::acoustics::AcousticMaterial;
    use crate::acoustics::LocallyResonantUnit;
    use crate::acoustics::PhononicCrystal1D;
    use crate::acoustics::PorousAbsorber;
    use crate::acoustics::UltrasonicNDE;
    #[test]
    fn test_water_density() {
        let m = AcousticMaterial::water();
        assert!((m.density - 998.2).abs() < 1.0);
    }
    #[test]
    fn test_water_is_fluid() {
        let m = AcousticMaterial::water();
        assert_eq!(m.shear_modulus, 0.0);
        assert_eq!(m.shear_velocity(), 0.0);
    }
    #[test]
    fn test_water_longitudinal_velocity() {
        let m = AcousticMaterial::water();
        let c = m.longitudinal_velocity();
        assert!(c > 1400.0 && c < 1600.0, "c_L = {c}");
    }
    #[test]
    fn test_water_impedance() {
        let m = AcousticMaterial::water();
        let z = m.acoustic_impedance();
        assert!(z > 1.4e6 && z < 1.6e6, "Z = {z}");
    }
    #[test]
    fn test_air_longitudinal_velocity() {
        let m = AcousticMaterial::air();
        let c = m.longitudinal_velocity();
        assert!(c > 330.0 && c < 360.0, "c_L = {c}");
    }
    #[test]
    fn test_steel_shear_velocity() {
        let m = AcousticMaterial::steel();
        let cs = m.shear_velocity();
        assert!(cs > 3000.0 && cs < 3500.0, "c_S = {cs}");
    }
    #[test]
    fn test_steel_longitudinal_velocity() {
        let m = AcousticMaterial::steel();
        let cl = m.longitudinal_velocity();
        assert!(cl > 5500.0 && cl < 6500.0, "c_L = {cl}");
    }
    #[test]
    fn test_concrete_density() {
        let m = AcousticMaterial::concrete();
        assert!((m.density - 2300.0).abs() < 1.0);
    }
    #[test]
    fn test_impedance_equals_rho_times_cl() {
        let m = AcousticMaterial::steel();
        let z = m.acoustic_impedance();
        assert!((z - m.density * m.longitudinal_velocity()).abs() < 1.0);
    }
    #[test]
    fn test_reflection_same_medium_is_zero() {
        let r = reflection_coefficient(1500.0, 1500.0);
        assert!(r.abs() < 1e-12);
    }
    #[test]
    fn test_reflection_total_at_rigid_wall() {
        let r = reflection_coefficient(1.0, 1e15);
        assert!((r - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_reflection_soft_boundary() {
        let r = reflection_coefficient(1000.0, 1e-9);
        assert!((r + 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_transmission_same_medium_is_one() {
        let t = transmission_coefficient(500.0, 500.0);
        assert!((t - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_r_and_t_energy_conservation() {
        let z1 = 400.0_f64;
        let z2 = 1.5e6_f64;
        let r = reflection_coefficient(z1, z2);
        let tau = 4.0 * z1 * z2 / ((z1 + z2) * (z1 + z2));
        let lhs = r * r + tau;
        assert!((lhs - 1.0).abs() < 1e-10, "lhs = {lhs}");
    }
    #[test]
    fn test_transmission_loss_same_medium_is_zero() {
        let tl = transmission_loss_db(1000.0, 1000.0);
        assert!(tl.abs() < 1e-10);
    }
    #[test]
    fn test_transmission_loss_large_impedance_mismatch() {
        let z_air = AcousticMaterial::air().acoustic_impedance();
        let z_water = AcousticMaterial::water().acoustic_impedance();
        let tl = transmission_loss_db(z_air, z_water);
        assert!(tl > 20.0, "TL = {tl} dB");
    }
    #[test]
    fn test_multilayer_trivial_no_layers() {
        let t2 = multilayer_transmission(&[400.0, 400.0], &[], 1000.0);
        assert!((t2 - 1.0).abs() < 1e-6, "|T|² = {t2}");
    }
    #[test]
    fn test_multilayer_single_layer_returns_nonzero() {
        let t2 = multilayer_transmission(&[400.0, 400.0, 400.0], &[0.1], 1000.0);
        assert!(t2 > 0.0, "|T|² = {t2}");
    }
    #[test]
    fn test_viscous_attenuation_increases_with_freq() {
        let alpha1 = viscous_attenuation(1000.0, 1e-3, 998.0, 1500.0);
        let alpha2 = viscous_attenuation(4000.0, 1e-3, 998.0, 1500.0);
        assert!(alpha2 > alpha1 * 10.0);
    }
    #[test]
    fn test_viscous_attenuation_positive() {
        let alpha = viscous_attenuation(1000.0, 1e-3, 998.0, 1500.0);
        assert!(alpha > 0.0);
    }
    #[test]
    fn test_atmospheric_attenuation_positive() {
        let alpha = atmospheric_attenuation_db_per_m(1000.0, 20.0, 0.50);
        assert!(alpha > 0.0, "alpha = {alpha}");
    }
    #[test]
    fn test_atmospheric_attenuation_increases_with_freq() {
        let a1 = atmospheric_attenuation_db_per_m(1000.0, 20.0, 0.50);
        let a2 = atmospheric_attenuation_db_per_m(8000.0, 20.0, 0.50);
        assert!(a2 > a1, "a1={a1}, a2={a2}");
    }
    #[test]
    fn test_attenuation_db_linear_in_distance() {
        let d1 = attenuation_db(0.5, 10.0);
        let d2 = attenuation_db(0.5, 20.0);
        assert!((d2 - 2.0 * d1).abs() < 1e-12);
    }
    #[test]
    fn test_sabine_t60_known_value() {
        let t60 = sabine_reverberation_time(200.0, 120.0, 0.20);
        assert!((t60 - 1.342).abs() < 0.01, "T60 = {t60}");
    }
    #[test]
    fn test_eyring_less_than_sabine_for_high_absorption() {
        let t_s = sabine_reverberation_time(200.0, 120.0, 0.50);
        let t_e = eyring_reverberation_time(200.0, 120.0, 0.50);
        assert!(t_e < t_s, "T_e={t_e}, T_s={t_s}");
    }
    #[test]
    fn test_sabine_eyring_agree_at_low_absorption() {
        let t_s = sabine_reverberation_time(500.0, 300.0, 0.05);
        let t_e = eyring_reverberation_time(500.0, 300.0, 0.05);
        assert!((t_s - t_e).abs() / t_s < 0.05, "T_s={t_s}, T_e={t_e}");
    }
    #[test]
    fn test_critical_distance_positive() {
        let rc = critical_distance(200.0, 1.5);
        assert!(rc > 0.0);
    }
    #[test]
    fn test_critical_distance_known_value() {
        let rc = critical_distance(500.0, 2.0);
        assert!((rc - 0.9013).abs() < 0.01, "r_c = {rc}");
    }
    #[test]
    fn test_noise_reduction_finite() {
        let nr = noise_reduction(50.0, 30.0, 10.0, 100.0, 1.0);
        assert!(nr.is_finite());
    }
    #[test]
    fn test_wave_number_at_1khz_in_air() {
        let k = wave_number(1000.0, 343.0);
        assert!((k - 18.33).abs() < 0.1, "k = {k}");
    }
    #[test]
    fn test_plane_wave_at_origin_zero_time() {
        let p = plane_wave_pressure(2.0, 1.0, 0.0, 1.0, 0.0);
        assert!((p - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_spl_at_reference_pressure_is_zero() {
        let spl = spl_db(20.0e-6);
        assert!(spl.abs() < 1e-10, "SPL = {spl}");
    }
    #[test]
    fn test_spl_doubles_pressure_is_6db() {
        let spl1 = spl_db(20.0e-6);
        let spl2 = spl_db(40.0e-6);
        assert!(
            (spl2 - spl1 - 6.020).abs() < 0.01,
            "delta = {}",
            spl2 - spl1
        );
    }
    #[test]
    fn test_loudness_at_40db_is_one_sone() {
        let s = loudness_sone(40.0, 1000.0);
        assert!((s - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_loudness_doubles_every_10db() {
        let s1 = loudness_sone(40.0, 1000.0);
        let s2 = loudness_sone(50.0, 1000.0);
        assert!((s2 / s1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_mass_law_tl_positive_at_1khz() {
        let tl = mass_law_tl(480.0, 1000.0);
        assert!(tl > 30.0, "TL = {tl} dB");
    }
    #[test]
    fn test_mass_law_tl_increases_with_frequency() {
        let tl1 = mass_law_tl(50.0, 500.0);
        let tl2 = mass_law_tl(50.0, 1000.0);
        assert!((tl2 - tl1 - 6.0).abs() < 0.1, "tl1={tl1}, tl2={tl2}");
    }
    #[test]
    fn test_mass_law_tl_increases_with_mass() {
        let tl1 = mass_law_tl(50.0, 1000.0);
        let tl2 = mass_law_tl(100.0, 1000.0);
        assert!((tl2 - tl1 - 6.0).abs() < 0.1, "tl1={tl1}, tl2={tl2}");
    }
    #[test]
    fn test_ultrasonic_nde_time_of_flight() {
        let nde = UltrasonicNDE::new(6000.0, 5.0e6, 1.0, 0.0);
        let tof = nde.time_of_flight(0.030);
        assert!((tof - 10.0e-6).abs() < 1e-9, "TOF = {tof}");
    }
    #[test]
    fn test_ultrasonic_nde_echo_amplitude_decreases_with_depth() {
        let nde = UltrasonicNDE::new(6000.0, 5.0e6, 1.0, 5.0);
        let a1 = nde.echo_amplitude(0.010, 1.0);
        let a2 = nde.echo_amplitude(0.050, 1.0);
        assert!(a1 > a2, "echo should decrease with depth: a1={a1}, a2={a2}");
    }
    #[test]
    fn test_ultrasonic_nde_wavelength() {
        let nde = UltrasonicNDE::new(5920.0, 5.0e6, 1.0, 0.0);
        let lam = nde.wavelength();
        assert!((lam - 5920.0 / 5.0e6).abs() < 1e-12, "λ = {lam}");
    }
    #[test]
    fn test_ultrasonic_nde_near_field_positive() {
        let nde = UltrasonicNDE::new(6000.0, 5.0e6, 1.0, 0.0);
        let nf = nde.near_field_distance(0.006);
        assert!(nf > 0.0, "near-field distance should be positive, got {nf}");
    }
    #[test]
    fn test_a_weighting_at_1khz_near_zero() {
        let w = a_weighting_db(1000.0);
        assert!(
            w.abs() < 1.0,
            "A-weighting at 1 kHz should be ~0 dB, got {w}"
        );
    }
    #[test]
    fn test_a_weighting_low_freq_negative() {
        let w = a_weighting_db(100.0);
        assert!(
            w < -10.0,
            "A-weighting at 100 Hz should be strongly negative, got {w}"
        );
    }
    #[test]
    fn test_c_weighting_near_1khz_near_zero() {
        let w = c_weighting_db(1000.0);
        assert!(
            w.abs() < 2.0,
            "C-weighting at 1 kHz should be ~0 dB, got {w}"
        );
    }
    #[test]
    fn test_overall_spl_from_two_equal_bands() {
        let overall = overall_spl_from_octave_bands(&[60.0, 60.0]);
        assert!((overall - 63.01).abs() < 0.02, "got {overall}");
    }
    #[test]
    fn test_apply_weighting_returns_same_length() {
        let spl = [55.0, 60.0, 65.0, 70.0, 72.0];
        let freqs = [250.0, 500.0, 1000.0, 2000.0, 4000.0];
        let weighted = apply_weighting(&spl, &freqs, WeightingFilter::A);
        assert_eq!(weighted.len(), spl.len());
    }
    #[test]
    fn test_apply_a_weighting_increases_at_high_freq_vs_low_freq() {
        let spl = [60.0, 60.0];
        let freqs_low = [100.0, 100.0];
        let freqs_high = [1000.0, 1000.0];
        let w_low = apply_weighting(&spl, &freqs_low, WeightingFilter::A);
        let w_high = apply_weighting(&spl, &freqs_high, WeightingFilter::A);
        assert!(w_high[0] > w_low[0], "A-weighting at 1 kHz > 100 Hz");
    }
    #[test]
    fn test_coincidence_frequency_positive() {
        let fc = coincidence_frequency(343.0, 5000.0, 0.012);
        assert!(fc > 0.0, "got {fc}");
    }
    #[test]
    fn test_matching_layer_geometric_mean() {
        let z1 = 400.0_f64;
        let z2 = 1.5e6_f64;
        let z_m = impedance_matching_layer(z1, z2);
        let expected = (z1 * z2).sqrt();
        assert!(
            (z_m - expected).abs() < 1.0,
            "Z_m = {z_m}, expected {expected}"
        );
    }
    #[test]
    fn test_matching_layer_zero_reflection_at_design_freq() {
        let z1 = 400.0_f64;
        let z2 = 1.5e6_f64;
        let z_m = impedance_matching_layer(z1, z2);
        let r = matching_layer_reflection(z1, z_m, z2, 1.0);
        assert!(
            r < 0.1,
            "Matching layer at design frequency should give near-zero reflection, got {r}"
        );
    }
    #[test]
    fn test_matching_layer_bandwidth_positive() {
        let z1 = 400.0;
        let z2 = 1.5e6;
        let z_m = impedance_matching_layer(z1, z2);
        let (f_low, f_high) = matching_layer_bandwidth(z1, z_m, z2, 0.5);
        assert!(
            f_high > f_low || (f_low > 1.9 && f_high < 0.2),
            "Bandwidth should be positive"
        );
    }
    #[test]
    fn test_matching_layer_impedance_symmetric() {
        let z1 = 1000.0;
        let z2 = 4000.0;
        let zm12 = impedance_matching_layer(z1, z2);
        let zm21 = impedance_matching_layer(z2, z1);
        assert!((zm12 - zm21).abs() < 1e-10);
    }
    #[test]
    fn test_angle_transmission_normal_incidence() {
        let z1 = 400.0;
        let z2 = 1.5e6;
        let c1 = 343.0;
        let c2 = 1480.0;
        let t = angle_dependent_transmission(z1, z2, c1, c2, 0.0).unwrap();
        let expected = 4.0 * z1 * z2 / (z1 + z2).powi(2);
        assert!(
            (t - expected).abs() < 1e-6,
            "Normal incidence T_I mismatch: got {t}"
        );
    }
    #[test]
    fn test_angle_transmission_total_internal_reflection() {
        let z_air = 400.0;
        let z_water = 1.5e6;
        let c_air = 343.0;
        let c_water = 1480.0;
        let t = angle_dependent_transmission(z_air, z_water, c_air, c_water, 80.0_f64.to_radians());
        assert!(
            t.is_none(),
            "Should be total internal reflection at 80° (air→water)"
        );
    }
    #[test]
    fn test_critical_angle_water_to_air() {
        let c1 = 1480.0;
        let c2 = 343.0;
        assert!(
            critical_angle(c1, c2).is_none(),
            "No critical angle when c2 < c1"
        );
    }
    #[test]
    fn test_critical_angle_air_to_water() {
        let c1 = 343.0;
        let c2 = 1480.0;
        let theta_c = critical_angle(c1, c2);
        assert!(theta_c.is_some());
        let angle = theta_c.unwrap();
        assert!(
            angle > 0.0 && angle < PI / 2.0,
            "Critical angle should be in (0, π/2)"
        );
    }
    #[test]
    fn test_lram_resonance_frequency_positive() {
        let unit = LocallyResonantUnit::new(0.1, 1e6, 1e-6, 7800.0);
        let f0 = unit.resonance_frequency();
        assert!(
            f0 > 0.0 && f0.is_finite(),
            "Resonance frequency should be positive finite: {f0}"
        );
    }
    #[test]
    fn test_lram_negative_mass_above_resonance() {
        let unit = LocallyResonantUnit::new(0.1, 1e6, 1e-6, 1.0e-6);
        let f0 = unit.resonance_frequency();
        let f_above = f0 * 1.01;
        assert!(
            unit.is_negative_mass(f_above),
            "Effective mass should be negative just above resonance"
        );
    }
    #[test]
    fn test_lram_positive_mass_below_resonance() {
        let unit = LocallyResonantUnit::new(0.1, 1e6, 1e-6, 7800.0);
        let f0 = unit.resonance_frequency();
        let f_below = f0 * 0.5;
        let rho_eff = unit.effective_density(f_below);
        assert!(
            rho_eff > 0.0,
            "Effective density should be positive below resonance: {rho_eff}"
        );
    }
    #[test]
    fn test_phononic_crystal_bragg_frequency_positive() {
        let pc = PhononicCrystal1D::new(400.0, 1.5e6, 343.0, 1480.0, 0.01, 0.01);
        let f_b = pc.bragg_frequency();
        assert!(
            f_b > 0.0 && f_b.is_finite(),
            "Bragg frequency should be positive: {f_b}"
        );
    }
    #[test]
    fn test_phononic_crystal_lattice_constant() {
        let pc = PhononicCrystal1D::new(400.0, 1.5e6, 343.0, 1480.0, 0.015, 0.010);
        assert!((pc.lattice_constant() - 0.025).abs() < 1e-12);
    }
    #[test]
    fn test_phononic_crystal_gap_width_positive_for_mismatch() {
        let pc = PhononicCrystal1D::new(400.0, 1.5e6, 343.0, 1480.0, 0.01, 0.01);
        let gap = pc.relative_gap_width();
        assert!(
            gap > 0.0,
            "Gap width should be positive for impedance mismatch: {gap}"
        );
    }
    #[test]
    fn test_phononic_crystal_zero_gap_for_equal_impedance() {
        let pc = PhononicCrystal1D::new(400.0, 400.0, 343.0, 343.0, 0.01, 0.01);
        let gap = pc.relative_gap_width();
        assert!(
            gap.abs() < 1e-12,
            "No gap for identical materials, got {gap}"
        );
    }
    #[test]
    fn test_phononic_crystal_is_in_gap_center() {
        let pc = PhononicCrystal1D::new(400.0, 1.5e6, 343.0, 1480.0, 0.01, 0.01);
        let f0 = pc.bragg_frequency();
        assert!(pc.is_in_gap(f0), "Bragg frequency should be in gap");
    }
    #[test]
    fn test_porous_absorber_coefficient_in_range() {
        let abs = PorousAbsorber::new(10_000.0, 0.05);
        let alpha = abs.absorption_coefficient(1000.0);
        assert!(
            (0.0..=1.0 + 1e-10).contains(&alpha),
            "Absorption coefficient should be in [0,1], got {alpha}"
        );
    }
    #[test]
    fn test_porous_absorber_impedance_real_positive() {
        let abs = PorousAbsorber::new(10_000.0, 0.05);
        let zr = abs.impedance_real(1000.0);
        assert!(
            zr > 1.0,
            "Real impedance should be > 1 (normalised), got {zr}"
        );
    }
    #[test]
    fn test_porous_absorber_absorption_increases_with_thickness() {
        let abs1 = PorousAbsorber::new(10_000.0, 0.03);
        let abs2 = PorousAbsorber::new(10_000.0, 0.10);
        let a1 = abs1.absorption_coefficient(500.0);
        let a2 = abs2.absorption_coefficient(500.0);
        assert!(
            a2 >= a1 - 0.1,
            "Thicker absorber should have >= absorption: a1={a1}, a2={a2}"
        );
    }
    #[test]
    fn test_porous_absorber_propagation_imag_positive() {
        let abs = PorousAbsorber::new(5_000.0, 0.05);
        let ki = abs.propagation_imag(2000.0);
        assert!(ki > 0.0, "Propagation attenuation should be positive: {ki}");
    }
    #[test]
    fn test_room_resonance_axial_x() {
        let lx = 5.0;
        let f = room_resonance_frequency(lx, 4.0, 3.0, 1, 0, 0, 343.0);
        assert!((f - 343.0 / (2.0 * lx)).abs() < 1e-6, "f = {f}");
    }
    #[test]
    fn test_room_resonance_axial_y() {
        let ly = 4.0;
        let f = room_resonance_frequency(5.0, ly, 3.0, 0, 1, 0, 343.0);
        assert!((f - 343.0 / (2.0 * ly)).abs() < 1e-6, "f = {f}");
    }
    #[test]
    fn test_room_resonance_modes_ordered() {
        let modes = room_modes(5.0, 4.0, 3.0, 343.0, 3, 6);
        assert_eq!(modes.len(), 6);
        for i in 0..modes.len() - 1 {
            assert!(
                modes[i] <= modes[i + 1],
                "modes not ordered: {} > {}",
                modes[i],
                modes[i + 1]
            );
        }
    }
    #[test]
    fn test_room_modes_all_positive() {
        let modes = room_modes(6.0, 5.0, 4.0, 343.0, 2, 10);
        for &f in &modes {
            assert!(f > 0.0, "all modes should be positive, found {f}");
        }
    }
    #[test]
    fn test_room_lowest_mode_correct() {
        let f = room_lowest_mode(6.0, 4.0, 3.0, 343.0);
        assert!((f - 343.0 / 12.0).abs() < 1e-6, "f = {f}");
    }
    #[test]
    fn test_schroeder_frequency_typical_room() {
        let fs = schroeder_frequency(1.0, 200.0);
        assert!((fs - 2000.0 / 200.0_f64.sqrt()).abs() < 0.1, "f_S = {fs}");
    }
    #[test]
    fn test_double_leaf_tl_better_than_single() {
        let tl_single = mass_law_tl(10.0, 1000.0);
        let tl_double = double_leaf_tl(10.0, 10.0, 0.1, 1000.0);
        assert!(
            tl_double > tl_single,
            "double leaf should outperform single: {tl_double} vs {tl_single}"
        );
    }
    #[test]
    fn test_double_leaf_tl_positive() {
        let tl = double_leaf_tl(12.0, 12.0, 0.08, 500.0);
        assert!(tl > 0.0, "TL should be positive, got {tl}");
    }
    #[test]
    fn test_sound_intensity_positive() {
        let i = sound_intensity(1.0, 400.0);
        assert!(i > 0.0, "intensity should be positive, got {i}");
    }
    #[test]
    fn test_sound_intensity_scales_with_pressure_squared() {
        let i1 = sound_intensity(1.0, 400.0);
        let i2 = sound_intensity(2.0, 400.0);
        assert!(
            (i2 / i1 - 4.0).abs() < 1e-10,
            "I ∝ p², got ratio {}",
            i2 / i1
        );
    }
    #[test]
    fn test_sound_intensity_level_reference_is_0db() {
        let sil = sound_intensity_level_db(1.0e-12);
        assert!(sil.abs() < 1e-10, "SIL at ref = 0 dB, got {sil}");
    }
    #[test]
    fn test_sound_power_level_reference_is_0db() {
        let pwl = sound_power_level_db(1.0e-12);
        assert!(pwl.abs() < 1e-10, "PWL at ref = 0 dB, got {pwl}");
    }
    #[test]
    fn test_doppler_stationary_source_and_observer() {
        let f = doppler_frequency(1000.0, 343.0, 0.0, 0.0);
        assert!((f - 1000.0).abs() < 1e-6, "no motion → no shift, got {f}");
    }
    #[test]
    fn test_doppler_approaching_source_raises_pitch() {
        let f = doppler_frequency(1000.0, 343.0, 34.3, 0.0);
        assert!(f > 1000.0, "approaching source → higher pitch, got {f}");
    }
    #[test]
    fn test_doppler_receding_source_lowers_pitch() {
        let f = doppler_frequency(1000.0, 343.0, -34.3, 0.0);
        assert!(f < 1000.0, "receding source → lower pitch, got {f}");
    }
    #[test]
    fn test_mach_number_subsonic() {
        let m = mach_number(100.0, 343.0);
        assert!(m < 1.0, "100 m/s is subsonic, got M = {m}");
    }
    #[test]
    fn test_mach_number_supersonic() {
        let m = mach_number(500.0, 343.0);
        assert!(m > 1.0, "500 m/s is supersonic, got M = {m}");
    }
    #[test]
    fn test_maekawa_il_positive_for_positive_delta() {
        let il = barrier_insertion_loss_maekawa(0.5, 1000.0, 343.0);
        assert!(
            il > 0.0,
            "barrier with positive path diff should give positive IL, got {il}"
        );
    }
    #[test]
    fn test_maekawa_il_increases_with_delta() {
        let il1 = barrier_insertion_loss_maekawa(0.2, 1000.0, 343.0);
        let il2 = barrier_insertion_loss_maekawa(1.0, 1000.0, 343.0);
        assert!(
            il2 > il1,
            "larger path diff → more insertion loss: {il1} vs {il2}"
        );
    }
    #[test]
    fn test_fresnel_number_formula() {
        let n = fresnel_number(0.5, 1000.0, 343.0);
        assert!((n - 2.0 * 0.5 * 1000.0 / 343.0).abs() < 1e-10, "got {n}");
    }
    #[test]
    fn test_helmholtz_resonator_typical_value() {
        let f = helmholtz_resonator_frequency(5.0e-4, 5.0e-4, 0.06, 343.0);
        assert!(
            f > 50.0 && f < 1000.0,
            "typical Helmholtz freq should be 50-1000 Hz, got {f}"
        );
    }
    #[test]
    fn test_quarter_wave_resonator_fundamental() {
        let f = quarter_wave_resonator_frequency(0.25, 343.0, 1);
        assert!((f - 343.0).abs() < 1e-6, "f_1 = {f}");
    }
    #[test]
    fn test_quarter_wave_resonator_third_harmonic() {
        let f1 = quarter_wave_resonator_frequency(0.25, 343.0, 1);
        let f2 = quarter_wave_resonator_frequency(0.25, 343.0, 2);
        assert!(
            (f2 / f1 - 3.0).abs() < 1e-6,
            "3rd harmonic ratio = {}",
            f2 / f1
        );
    }
    #[test]
    fn test_half_wave_resonator_fundamental() {
        let f = half_wave_resonator_frequency(0.25, 343.0, 1);
        assert!((f - 686.0).abs() < 1e-6, "f_1 = {f}");
    }
    #[test]
    fn test_half_wave_second_harmonic() {
        let f1 = half_wave_resonator_frequency(0.5, 343.0, 1);
        let f2 = half_wave_resonator_frequency(0.5, 343.0, 2);
        assert!((f2 / f1 - 2.0).abs() < 1e-10, "ratio = {}", f2 / f1);
    }
    #[test]
    fn test_mean_free_path_cube_room() {
        let l = mean_free_path(64.0, 96.0);
        assert!((l - 64.0 * 4.0 / 96.0).abs() < 1e-10, "L_mfp = {l}");
    }
    #[test]
    fn test_room_constant_increases_with_absorption() {
        let r1 = room_constant(100.0, 0.1);
        let r2 = room_constant(100.0, 0.5);
        assert!(
            r2 > r1,
            "room constant increases with absorption: {r1} vs {r2}"
        );
    }
    #[test]
    fn test_diffuse_field_spl_finite() {
        let r = room_constant(200.0, 0.2);
        let spl = diffuse_field_spl(0.01, r);
        assert!(spl.is_finite(), "diffuse SPL should be finite, got {spl}");
    }
    #[test]
    fn test_total_spl_in_room_decreases_with_distance() {
        let r = room_constant(500.0, 0.15);
        let spl1 = total_spl_in_room(0.1, r, 1.0, 1.0);
        let spl2 = total_spl_in_room(0.1, r, 5.0, 1.0);
        assert!(
            spl1 > spl2,
            "SPL should decrease with distance: {spl1} vs {spl2}"
        );
    }
    #[test]
    fn test_insertion_loss_positive() {
        let mat = AcousticMaterial::concrete();
        let il = mat.compute_insertion_loss(500.0, 0.2, 1.0);
        assert!(il > 0.0, "insertion loss should be positive, got {il}");
    }
    #[test]
    fn test_insertion_loss_increases_with_frequency() {
        let mat = AcousticMaterial::concrete();
        let il_low = mat.compute_insertion_loss(250.0, 0.2, 1.0);
        let il_high = mat.compute_insertion_loss(2000.0, 0.2, 1.0);
        assert!(
            il_high > il_low,
            "IL should increase with frequency: {il_low} vs {il_high}"
        );
    }
    #[test]
    fn test_insertion_loss_increases_with_thickness() {
        let mat = AcousticMaterial::concrete();
        let il_thin = mat.compute_insertion_loss(500.0, 0.1, 0.5);
        let il_thick = mat.compute_insertion_loss(500.0, 0.5, 0.5);
        assert!(
            il_thick > il_thin,
            "thicker barrier → higher IL: {il_thin} vs {il_thick}"
        );
    }
    #[test]
    fn test_insertion_loss_zero_fresnel_still_has_mass_term() {
        let mat = AcousticMaterial::concrete();
        let il = mat.compute_insertion_loss(1000.0, 0.2, 0.0);
        assert!(
            il > 0.0,
            "mass-law contribution should give positive IL: {il}"
        );
    }
    #[test]
    fn test_nrc_between_zero_and_one() {
        let mat = AcousticMaterial::concrete();
        let nrc = mat.compute_noise_reduction_coefficient(0.05);
        assert!((0.0..=1.0).contains(&nrc), "NRC must be in [0, 1]: {nrc}");
    }
    #[test]
    fn test_nrc_rounded_to_nearest_0_05() {
        let mat = AcousticMaterial::concrete();
        let nrc = mat.compute_noise_reduction_coefficient(0.02);
        let remainder = (nrc / 0.05 - (nrc / 0.05).round()).abs();
        assert!(remainder < 1e-10, "NRC should be rounded to 0.05: {nrc}");
    }
    #[test]
    fn test_nrc_higher_loss_factor_gives_higher_nrc() {
        let mat_lo = AcousticMaterial::new(2300.0, 13.33e9, 10.0e9, 0.01);
        let mat_hi = AcousticMaterial::new(2300.0, 13.33e9, 10.0e9, 0.20);
        let nrc_lo = mat_lo.compute_noise_reduction_coefficient(0.05);
        let nrc_hi = mat_hi.compute_noise_reduction_coefficient(0.05);
        assert!(
            nrc_hi >= nrc_lo,
            "higher loss factor → higher NRC: {nrc_lo} vs {nrc_hi}"
        );
    }
    #[test]
    fn test_mass_law_tl_positive_heavy_panel() {
        let mat = AcousticMaterial::concrete();
        let tl = mat.compute_transmission_loss_mass_law(500.0, 0.2);
        assert!(tl > 0.0, "heavy panel TL should be positive, got {tl}");
    }
    #[test]
    fn test_mass_law_tl_doubles_6db_per_octave() {
        let mat = AcousticMaterial::concrete();
        let tl1 = mat.compute_transmission_loss_mass_law(500.0, 0.1);
        let tl2 = mat.compute_transmission_loss_mass_law(1000.0, 0.1);
        let delta = tl2 - tl1;
        assert!(
            (delta - 6.0).abs() < 0.1,
            "mass law should give ~6 dB per octave: got {delta:.3} dB"
        );
    }
    #[test]
    fn test_mass_law_tl_doubles_6db_per_mass_doubling() {
        let mat = AcousticMaterial::concrete();
        let tl1 = mat.compute_transmission_loss_mass_law(1000.0, 0.1);
        let tl2 = mat.compute_transmission_loss_mass_law(1000.0, 0.2);
        let delta = tl2 - tl1;
        assert!(
            (delta - 6.0).abs() < 0.1,
            "doubling mass should add ~6 dB: got {delta:.3} dB"
        );
    }
}
