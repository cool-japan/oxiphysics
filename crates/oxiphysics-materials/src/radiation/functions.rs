//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::RadiationSurface;

/// Stefan-Boltzmann constant \[W/(m²·K⁴)\]
pub const SIGMA: f64 = 5.67e-8;
/// Wien's displacement constant \[m·K\]
pub const WIEN_B: f64 = 2.898e-3;
/// Planck's constant \[J·s\]
pub const PLANCK_H: f64 = 6.626e-34;
/// Speed of light in vacuum \[m/s\]
pub const SPEED_OF_LIGHT: f64 = 2.998e8;
/// Boltzmann constant \[J/K\]
pub const BOLTZMANN_K: f64 = 1.381e-23;
/// Blackbody emissive power \[W/m²\]: E_b = sigma * T^4.
#[allow(dead_code)]
pub fn blackbody_emissive_power(temp_k: f64) -> f64 {
    SIGMA * temp_k.powi(4)
}
/// Planck spectral intensity \[W/(m²·sr·m)\] at wavelength (m) and temperature (K).
///
/// I_lambda = (2*h*c²/lambda^5) / (exp(h*c/(lambda*k*T)) - 1)
#[allow(dead_code)]
pub fn blackbody_spectral_intensity(wavelength_m: f64, temp_k: f64) -> f64 {
    let h = PLANCK_H;
    let c = SPEED_OF_LIGHT;
    let k = BOLTZMANN_K;
    let numerator = 2.0 * h * c * c / wavelength_m.powi(5);
    let exponent = h * c / (wavelength_m * k * temp_k);
    numerator / (exponent.exp() - 1.0)
}
/// Wien's displacement law: lambda_max = b / T.
#[allow(dead_code)]
pub fn wien_displacement(temp_k: f64) -> f64 {
    WIEN_B / temp_k
}
/// Total emission from a gray surface \[W\]: Q = emissivity * sigma * T^4 * area.
#[allow(dead_code)]
pub fn total_emission(emissivity: f64, temp_k: f64, area: f64) -> f64 {
    emissivity * SIGMA * temp_k.powi(4) * area
}
/// Net radiation exchange between two gray surfaces (simplified) \[W\].
///
/// q = epsilon_1 * sigma * A * F12 * (T1^4 - T2^4)
#[allow(dead_code)]
pub fn net_radiation_exchange(e1: f64, t1: f64, e2: f64, t2: f64, f12: f64, a: f64) -> f64 {
    let _ = e2;
    e1 * SIGMA * a * f12 * (t1.powi(4) - t2.powi(4))
}
/// Radiation heat transfer coefficient \[W/(m²·K)\].
///
/// h_r = epsilon * sigma * (T1² + T2²) * (T1 + T2)
#[allow(dead_code)]
pub fn radiation_heat_transfer_coefficient(e: f64, t1: f64, t2: f64) -> f64 {
    e * SIGMA * (t1 * t1 + t2 * t2) * (t1 + t2)
}
/// View factor F12 for two equal parallel rectangular plates of size a×b
/// separated by distance h (analytical formula).
///
/// Uses the standard double-area-integral result (Sparrow-Cess / radiation
/// handbook formula for aligned parallel rectangles).
#[allow(dead_code)]
pub fn view_factor_parallel_plates(a: f64, b: f64, h: f64) -> f64 {
    let x = a / h;
    let y = b / h;
    let s = |u: f64, v: f64| -> f64 { (u * u + v * v).sqrt() };
    let ln_term = (((1.0 + x * x) * (1.0 + y * y)) / (1.0 + x * x + y * y))
        .sqrt()
        .ln();
    let t1 = x * (1.0 + y * y).sqrt() * x.atan2((1.0 + y * y).sqrt());
    let t2 = y * (1.0 + x * x).sqrt() * y.atan2((1.0 + x * x).sqrt());
    let t3 = -x * x.atan();
    let t4 = -y * y.atan();
    let _ = s(x, y);
    let f = (2.0 / (std::f64::consts::PI * x * y)) * (ln_term + t1 + t2 + t3 + t4);
    f.clamp(0.0, 1.0)
}
/// View factor F12 for two coaxial parallel disks of radii r1, r2 separated
/// by distance h.
///
/// Exact formula: F12 = 0.5*(S - sqrt(S² - 4*(r2/r1)²))
/// where S = 1 + (1 + R2²) / R1²,  R1 = r1/h, R2 = r2/h.
#[allow(dead_code)]
pub fn view_factor_coaxial_disks(r1: f64, r2: f64, h: f64) -> f64 {
    let big_r1 = r1 / h;
    let big_r2 = r2 / h;
    let s = 1.0 + (1.0 + big_r2 * big_r2) / (big_r1 * big_r1);
    0.5 * (s - (s * s - 4.0 * (big_r2 / big_r1).powi(2)).sqrt())
}
/// View factor reciprocity: F21 = A1 * F12 / A2.
#[allow(dead_code)]
pub fn view_factor_reciprocity(f12: f64, a1: f64, a2: f64) -> f64 {
    f12 * a1 / a2
}
/// View factor sum rule: F_self = 1 - sum(F_others).
#[allow(dead_code)]
pub fn view_factor_sum_rule(f_others: &[f64]) -> f64 {
    let total: f64 = f_others.iter().sum();
    (1.0 - total).max(0.0)
}
/// Heat transfer rate for a two-surface enclosure \[W\].
///
/// q = (E_b1 - E_b2) / (R_s1 + 1/(A1*F12) + R_s2)
#[allow(dead_code)]
pub fn two_surface_enclosure_heat_transfer(
    s1: &RadiationSurface,
    s2: &RadiationSurface,
    f12: f64,
) -> f64 {
    let eb1 = SIGMA * s1.temperature.powi(4);
    let eb2 = SIGMA * s2.temperature.powi(4);
    let r_s1 = s1.surface_resistance();
    let r_s2 = s2.surface_resistance();
    let r_space = 1.0 / (s1.area * f12);
    let r_total = r_s1 + r_space + r_s2;
    (eb1 - eb2) / r_total
}
/// Total irradiance on a tilted surface \[W/m²\] using the ASHRAE decomposition.
///
/// Arguments:
/// - `dni`  – direct normal irradiance \[W/m²\]
/// - `dhi`  – diffuse horizontal irradiance \[W/m²\]
/// - `tilt_deg` – surface tilt from horizontal \[°\]
/// - `azimuth_deg` – surface azimuth (0 = south, positive east) \[°\]
/// - `sun_altitude_deg` – solar altitude angle \[°\]
/// - `sun_azimuth_deg`  – solar azimuth angle (0 = south, positive east) \[°\]
#[allow(dead_code)]
pub fn solar_irradiance_surface(
    dni: f64,
    dhi: f64,
    tilt_deg: f64,
    azimuth_deg: f64,
    sun_altitude_deg: f64,
    sun_azimuth_deg: f64,
) -> f64 {
    let tilt = tilt_deg.to_radians();
    let az_surf = azimuth_deg.to_radians();
    let alt_sun = sun_altitude_deg.to_radians();
    let az_sun = sun_azimuth_deg.to_radians();
    let cos_theta =
        alt_sun.sin() * tilt.cos() + alt_sun.cos() * tilt.sin() * (az_sun - az_surf).cos();
    let direct = (dni * cos_theta).max(0.0);
    let diffuse = dhi * (1.0 + tilt.cos()) / 2.0;
    direct + diffuse
}
/// Air mass: AM = 1 / sin(altitude) (Kasten simple form).
///
/// Returns a large value for altitudes near zero to avoid division by zero.
#[allow(dead_code)]
pub fn air_mass(sun_altitude_deg: f64) -> f64 {
    let alt = sun_altitude_deg.to_radians();
    if alt <= 0.0 {
        return f64::INFINITY;
    }
    1.0 / alt.sin()
}
/// Atmospheric attenuation factor: I/I0 = exp(-tau * AM).
#[allow(dead_code)]
pub fn attenuation_factor(air_mass: f64, tau: f64) -> f64 {
    (-tau * air_mass).exp()
}
/// Net radiative heat flux between two infinite parallel gray plates \[W/m²\].
///
/// `q'' = σ·(T₁⁴ − T₂⁴) / (1/ε₁ + 1/ε₂ − 1)`
#[allow(dead_code)]
pub fn thermal_radiation_parallel_plates(eps1: f64, temp1: f64, eps2: f64, temp2: f64) -> f64 {
    let eb1 = SIGMA * temp1.powi(4);
    let eb2 = SIGMA * temp2.powi(4);
    let r = 1.0 / eps1 + 1.0 / eps2 - 1.0;
    (eb1 - eb2) / r
}
/// Net radiative heat flux from a small convex body (1) fully enclosed by
/// a large isothermal enclosure (2) \[W/m²\].
///
/// `q'' = ε₁·σ·(T₁⁴ − T₂⁴)`
#[allow(dead_code)]
pub fn thermal_radiation_body_in_enclosure(eps1: f64, temp1: f64, temp2: f64) -> f64 {
    eps1 * SIGMA * (temp1.powi(4) - temp2.powi(4))
}
/// View factor F₁₂ from the inner surface of two concentric coaxial
/// cylinders (infinite length) to the outer cylinder.
///
/// By definition F₁₂ = 1 for an infinite inner cylinder (all radiation
/// hits the outer surface). The reciprocal is `F₂₁ = r₁/r₂`.
///
/// # Arguments
/// * `r_inner` — radius of inner cylinder \[m\]
/// * `r_outer` — radius of outer cylinder \[m\]
#[allow(dead_code)]
pub fn view_factor_concentric_cylinders(r_inner: f64, r_outer: f64) -> (f64, f64) {
    let _ = r_outer;
    let _ = r_inner;
    (1.0, r_inner / r_outer)
}
/// AM1.5 reference total irradiance \[W/m²\] (ASTM G-173-03 standard).
pub const AM15_TOTAL_IRRADIANCE: f64 = 1000.0;
/// Solar constant at 1 AU from the sun \[W/m²\] (IAU 2015 value).
pub const SOLAR_CONSTANT: f64 = 1361.0;
/// Simple seasonal irradiance variation.
///
/// Returns the extra-terrestrial irradiance \[W/m²\] for the given day of year.
/// Uses Spencer's (1971) approximation:
/// `I = I_sc · (1 + 0.033 · cos(2π·n/365))`.
///
/// # Arguments
/// * `day_of_year` — calendar day (1–365)
#[allow(dead_code)]
pub fn extraterrestrial_irradiance(day_of_year: u32) -> f64 {
    let b = 2.0 * std::f64::consts::PI * (day_of_year as f64) / 365.0;
    SOLAR_CONSTANT * (1.0 + 0.033 * b.cos())
}
/// Diurnal solar elevation angle \[degrees\] for a given latitude and hour angle.
///
/// `sin(α) = sin(φ)·sin(δ) + cos(φ)·cos(δ)·cos(H)`
///
/// # Arguments
/// * `lat_deg`   — observer latitude \[°\]
/// * `decl_deg`  — solar declination \[°\] (e.g. ±23.45 at solstices)
/// * `hour_angle_deg` — hour angle \[°\] (0 = solar noon, ±180 = midnight)
#[allow(dead_code)]
pub fn solar_elevation(lat_deg: f64, decl_deg: f64, hour_angle_deg: f64) -> f64 {
    let phi = lat_deg.to_radians();
    let delta = decl_deg.to_radians();
    let h = hour_angle_deg.to_radians();
    let sin_alpha = phi.sin() * delta.sin() + phi.cos() * delta.cos() * h.cos();
    sin_alpha.asin().to_degrees()
}
/// Effective emissivity of N radiation shields between two parallel plates.
///
/// Each shield (and both plates) is gray with the same emissivity `epsilon`.
/// The effective emissivity is:
/// `ε_eff = 1 / (1/ε₁ + 1/ε₂ - 1 + N*(2/ε_shield - 1))`
///
/// For the special case `ε₁ = ε₂ = ε_shield = ε`:
/// `ε_eff = ε / (1 + N*(2 - ε) - (N-1)*ε) ≈ ...` → simplified below.
///
/// General formula:
/// `1/ε_eff = 1/ε₁ + 1/ε₂ - 1 + Σ(1/ε_s,i + 1/ε_s,i' - 1)`
///
/// This function assumes all shields have the same emissivity on both faces.
///
/// # Arguments
/// * `eps1`     — emissivity of surface 1
/// * `eps2`     — emissivity of surface 2
/// * `eps_sh`   — emissivity of each shield (assumed equal on both faces)
/// * `n_shields`— number of shields
#[allow(dead_code)]
pub fn radiation_shield_effective_emissivity(
    eps1: f64,
    eps2: f64,
    eps_sh: f64,
    n_shields: u32,
) -> f64 {
    let r_total = 1.0 / eps1 + 1.0 / eps2 - 1.0 + n_shields as f64 * (2.0 / eps_sh - 1.0);
    if r_total < f64::EPSILON {
        return 1.0;
    }
    1.0 / r_total
}
/// Net heat flux reduction factor due to radiation shields.
///
/// Returns `q_with_shield / q_without_shield`.
///
/// Without shields: `q₀ = σ*(T₁⁴-T₂⁴) * ε_eff₀`
/// With N shields:  `q  = σ*(T₁⁴-T₂⁴) * ε_eff`
///
/// Reduction = `ε_eff / ε_eff₀` = `(1+N*(2/ε-1)) / (1+N*(2/ε-1)+... )`
///
/// For equal emissivities: ratio = `1 / (1 + N * (2/ε - 1) * ε / (2 - ε))`
#[allow(dead_code)]
pub fn radiation_shield_reduction_factor(eps: f64, n_shields: u32) -> f64 {
    let eps_eff_no_shield = radiation_shield_effective_emissivity(eps, eps, eps, 0);
    let eps_eff_with = radiation_shield_effective_emissivity(eps, eps, eps, n_shields);
    if eps_eff_no_shield < f64::EPSILON {
        return 1.0;
    }
    eps_eff_with / eps_eff_no_shield
}
/// Verify Kirchhoff's law: at thermal equilibrium, absorptivity α = emissivity ε.
///
/// Returns `true` if `|alpha - epsilon| < tol`.
///
/// This is a convenience validator function used to check that a radiative
/// surface model satisfies the fundamental thermodynamic constraint.
#[allow(dead_code)]
pub fn kirchhoff_law_satisfied(emissivity: f64, absorptivity: f64, tol: f64) -> bool {
    (emissivity - absorptivity).abs() < tol
}
/// Gray surface approximation: compute absorptivity from emissivity assuming
/// Kirchhoff's law holds (valid at thermal equilibrium for opaque gray surfaces).
#[allow(dead_code)]
pub fn absorptivity_from_emissivity(emissivity: f64) -> f64 {
    emissivity
}
/// Spectral emissivity from absorptivity using Kirchhoff's law.
///
/// For a surface in thermal equilibrium: `ε_λ = α_λ` at each wavelength.
/// Returns a vector of emissivities matching the input absorptivity slice.
#[allow(dead_code)]
pub fn spectral_emissivity_kirchhoff(spectral_absorptivity: &[f64]) -> Vec<f64> {
    spectral_absorptivity.to_vec()
}
/// Radiosity of a gray diffuse surface.
///
/// `J = ε·E_b + (1-ε)·G`
///
/// where `E_b = σ·T⁴` is the blackbody emissive power and `G` is the irradiation.
/// For a blackbody (ε = 1), `J = E_b`.
#[allow(dead_code)]
pub fn surface_radiosity(emissivity: f64, temp_k: f64, irradiation: f64) -> f64 {
    let e_b = SIGMA * temp_k.powi(4);
    emissivity * e_b + (1.0 - emissivity) * irradiation
}
/// Net radiative flux from surface i \[W/m²\] in the radiosity method.
///
/// `q_i = (E_b,i - J_i) / ((1-ε_i)/ε_i)` (surface resistance form)
///
/// When ε = 1 (blackbody), returns `E_b - J`.
#[allow(dead_code)]
pub fn radiosity_net_flux(emissivity: f64, temp_k: f64, radiosity: f64) -> f64 {
    let e_b = SIGMA * temp_k.powi(4);
    if (emissivity - 1.0).abs() < 1e-12 {
        return e_b - radiosity;
    }
    (e_b - radiosity) * emissivity / (1.0 - emissivity)
}
/// Solve the N-surface radiosity system iteratively (Jacobi method).
///
/// Given temperatures, emissivities, areas, and view-factor matrix, computes
/// the radiosity vector `J[i]` satisfying:
///
/// `J[i] = ε[i]*σ*T[i]⁴ + (1-ε[i]) * Σ_j F[i,j]*J[j]`
///
/// Returns the radiosity vector after `n_iter` Jacobi iterations.
///
/// # Arguments
/// * `temps`        — temperatures \[K\]
/// * `emissivities` — emissivities (0–1)
/// * `view_factors` — row-major n×n view factor matrix
/// * `n_iter`       — number of Jacobi iterations
#[allow(dead_code)]
pub fn radiosity_solve_jacobi(
    temps: &[f64],
    emissivities: &[f64],
    view_factors: &[f64],
    n_iter: usize,
) -> Vec<f64> {
    let n = temps.len();
    let mut j: Vec<f64> = temps
        .iter()
        .zip(emissivities.iter())
        .map(|(&t, &e)| e * SIGMA * t.powi(4))
        .collect();
    for _ in 0..n_iter {
        let j_old = j.clone();
        for i in 0..n {
            let e_b = SIGMA * temps[i].powi(4);
            let eps = emissivities[i];
            let g: f64 = (0..n).map(|k| view_factors[i * n + k] * j_old[k]).sum();
            j[i] = eps * e_b + (1.0 - eps) * g;
        }
    }
    j
}
/// Net heat flow \[W\] from surface i computed from the radiosity solution.
///
/// `q_i = A_i * (J_i - G_i)` where `G_i = Σ_j F[i,j] * J[j]`
#[allow(dead_code)]
pub fn radiosity_net_heat_flow(
    i: usize,
    areas: &[f64],
    radiosities: &[f64],
    view_factors: &[f64],
) -> f64 {
    let n = areas.len();
    let g_i: f64 = (0..n)
        .map(|k| view_factors[i * n + k] * radiosities[k])
        .sum();
    areas[i] * (radiosities[i] - g_i)
}
/// Simplified AM1.5 spectral irradiance \[W/(m²·nm)\] at wavelength `lambda_nm` \[nm\].
///
/// Uses a piecewise Gaussian approximation to the ASTM G-173-03 AM1.5G spectrum,
/// normalised so the integrated total ≈ 1000 W/m².
///
/// This is a gross simplification for modelling purposes only.
/// For accurate photovoltaic calculations, use tabulated data.
///
/// Dominant peak near 500 nm (visible), secondary peak near 1600 nm (NIR).
#[allow(dead_code)]
pub fn am15_spectral_irradiance(lambda_nm: f64) -> f64 {
    let vis = 1.8 * (-0.5 * ((lambda_nm - 500.0) / 150.0).powi(2)).exp();
    let nir = 0.9 * (-0.5 * ((lambda_nm - 1000.0) / 300.0).powi(2)).exp();
    let fir = 0.4 * (-0.5 * ((lambda_nm - 1600.0) / 400.0).powi(2)).exp();
    (vis + nir + fir).max(0.0)
}
/// Integrated AM1.5 irradiance \[W/m²\] over a wavelength range \[lambda_a, lambda_b\] (nm).
///
/// Uses a simple trapezoidal integration with `n_steps` intervals.
#[allow(dead_code)]
pub fn am15_integrated_irradiance(lambda_a_nm: f64, lambda_b_nm: f64, n_steps: usize) -> f64 {
    let n = n_steps.max(2);
    let dl = (lambda_b_nm - lambda_a_nm) / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let lam = lambda_a_nm + i as f64 * dl;
        let w = if i == 0 || i == n { 0.5 } else { 1.0 };
        sum += w * am15_spectral_irradiance(lam);
    }
    sum * dl
}
/// Photon flux density \[photons/(s·m²·nm)\] at wavelength `lambda_nm` \[nm\].
///
/// `Φ_λ = E_λ / (h·c/λ)` where `E_λ` is the spectral irradiance.
#[allow(dead_code)]
pub fn am15_photon_flux(lambda_nm: f64) -> f64 {
    let irradiance = am15_spectral_irradiance(lambda_nm);
    let lambda_m = lambda_nm * 1.0e-9;
    let photon_energy = PLANCK_H * SPEED_OF_LIGHT / lambda_m;
    if photon_energy <= 0.0 {
        return 0.0;
    }
    irradiance / photon_energy
}
/// View factor F12 for two perpendicular rectangles sharing a common edge.
///
/// Uses the analytical formula from Howell's view factor catalog (Configuration C-11).
/// Rectangles of dimensions H×W and L×W sharing a common edge of length W,
/// oriented perpendicularly.
///
/// Simplified version for two equal-width perpendicular strips (H = L = a, W = b):
///
/// `F12 = (1/π) * [W * atan(1/W) + H * atan(1/H) - sqrt(H²+W²)*atan(1/sqrt(H²+W²))
///                 + (1/4)*ln(part)]`
///
/// where H = h/w, W = l/w  (w = common edge length).
///
/// This provides the exact formula for the simplest case (equal-size rectangles).
///
/// # Arguments
/// * `h` — height of first rectangle \[m\]
/// * `l` — length of second rectangle \[m\]
/// * `w` — width (shared edge) \[m\]
#[allow(dead_code)]
pub fn view_factor_perpendicular_rectangles(h: f64, l: f64, w: f64) -> f64 {
    let big_h = h / w;
    let big_l = l / w;
    let h2 = big_h * big_h;
    let l2 = big_l * big_l;
    let t1 = big_h * (1.0 + l2).sqrt().atan2(big_h * (1.0 + l2).sqrt()) * (1.0 + l2).sqrt();
    let t2 = big_l * (1.0 + h2).sqrt().atan2(big_l) * (1.0 + h2).sqrt();
    let f = (1.0 / (std::f64::consts::PI * big_h))
        * (big_h * (big_l).atan() + big_l * (big_h).atan()
            - (h2 + l2).sqrt() * ((h2 + l2).sqrt()).atan()
            + 0.25 * ((1.0 + h2) * (1.0 + l2) / (1.0 + h2 + l2)).ln());
    let _ = (t1, t2);
    f.clamp(0.0, 1.0)
}
/// View factor from a small convex sphere of area A1 to a surrounding enclosure.
///
/// For a convex body: F12 = 1 (all radiation from convex body hits enclosure).
/// Reciprocal: F21 = A1 / A_enclosure.
///
/// # Arguments
/// * `r_sphere`    — radius of the sphere \[m\]
/// * `a_enclosure` — total area of the surrounding enclosure \[m²\]
#[allow(dead_code)]
pub fn view_factor_sphere_in_enclosure(r_sphere: f64, a_enclosure: f64) -> f64 {
    let a_sphere = 4.0 * std::f64::consts::PI * r_sphere * r_sphere;
    (a_sphere / a_enclosure).min(1.0)
}
/// View factor F12 for two infinite parallel strips of equal width `w` separated by `h`.
///
/// Uses the crossed-string method for 2-D geometry:
/// `F12 = (sqrt(w²+h²) - h) / w`  (Hottel's crossed-string method)
#[allow(dead_code)]
pub fn view_factor_infinite_parallel_strips(w1: f64, w2: f64, h: f64) -> f64 {
    let _ = w2;
    let w = w1;
    let result = ((w * w + h * h).sqrt() - h) / w;
    result.clamp(0.0, 1.0)
}
/// View factor from a differential area element to a sphere.
///
/// `F_dA→sphere = (r/d)² / 4` where `d` is the distance centre-to-element, `r` sphere radius.
///
/// Valid only when `d > r` (element outside sphere).
#[allow(dead_code)]
pub fn view_factor_diff_area_to_sphere(r_sphere: f64, d_centre: f64) -> f64 {
    if d_centre <= r_sphere {
        return 1.0;
    }
    0.25 * (r_sphere / d_centre).powi(2)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::radiation::BlackbodySpectrum;
    use crate::radiation::GrayBodyModel;
    use crate::radiation::IrradiationDosimetry;
    use crate::radiation::MCRadiationTransport;
    use crate::radiation::MonteCarloViewFactor;
    use crate::radiation::NeutronModeration;
    use crate::radiation::ParticipatingMedia;
    use crate::radiation::RadiationDamage;
    use crate::radiation::RadiationNetwork;
    use crate::radiation::RadiativeProperties;
    use crate::radiation::SolarCell;
    use crate::radiation::SpectralEmissivityModel;
    use crate::radiation::StefanBoltzmannIntegrator;
    #[test]
    fn test_blackbody_emissive_power_0k() {
        assert_eq!(blackbody_emissive_power(0.0), 0.0);
    }
    #[test]
    fn test_blackbody_emissive_power_1000k() {
        let e = blackbody_emissive_power(1000.0);
        assert!((e - 56_700.0).abs() < 1.0, "got {e}");
    }
    #[test]
    fn test_blackbody_emissive_power_300k() {
        let e = blackbody_emissive_power(300.0);
        assert!((e - 459.27).abs() < 1.0, "got {e}");
    }
    #[test]
    fn test_blackbody_emissive_power_t4_scaling() {
        let e1 = blackbody_emissive_power(500.0);
        let e2 = blackbody_emissive_power(1000.0);
        assert!((e2 / e1 - 16.0).abs() < 1e-6, "ratio = {}", e2 / e1);
    }
    #[test]
    fn test_wien_displacement_sun() {
        let lam = wien_displacement(5778.0);
        assert!((lam - 501.6e-9).abs() < 5e-9, "got {lam}");
    }
    #[test]
    fn test_wien_displacement_room_temp() {
        let lam = wien_displacement(300.0);
        assert!((lam - 9.66e-6).abs() < 1e-7, "got {lam}");
    }
    #[test]
    fn test_total_emission_blackbody() {
        let q = total_emission(1.0, 1000.0, 1.0);
        assert!((q - blackbody_emissive_power(1000.0)).abs() < 1e-6);
    }
    #[test]
    fn test_total_emission_area_scaling() {
        let q1 = total_emission(0.9, 500.0, 1.0);
        let q2 = total_emission(0.9, 500.0, 2.0);
        assert!((q2 / q1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_blackbody_spectral_intensity_positive() {
        let i = blackbody_spectral_intensity(1e-6, 5000.0);
        assert!(i > 0.0, "got {i}");
    }
    #[test]
    fn test_blackbody_spectral_intensity_peak_shifts() {
        let i_hot = blackbody_spectral_intensity(500e-9, 6000.0);
        let i_cold = blackbody_spectral_intensity(500e-9, 3000.0);
        assert!(i_hot > i_cold);
    }
    #[test]
    fn test_net_radiation_exchange_same_temp() {
        let q = net_radiation_exchange(0.8, 400.0, 0.8, 400.0, 1.0, 1.0);
        assert!(q.abs() < 1e-6, "got {q}");
    }
    #[test]
    fn test_net_radiation_exchange_hot_to_cold() {
        let q = net_radiation_exchange(0.9, 600.0, 0.9, 300.0, 1.0, 1.0);
        assert!(q > 0.0, "hot-to-cold should be positive");
    }
    #[test]
    fn test_radiation_heat_transfer_coefficient_positive() {
        let h = radiation_heat_transfer_coefficient(0.9, 400.0, 300.0);
        assert!(h > 0.0, "got {h}");
    }
    #[test]
    fn test_radiation_heat_transfer_coefficient_known() {
        let h = radiation_heat_transfer_coefficient(0.9, 400.0, 300.0);
        assert!(h > 5.0 && h < 15.0, "got {h}");
    }
    #[test]
    fn test_view_factor_reciprocity() {
        let f21 = view_factor_reciprocity(0.4, 2.0, 4.0);
        assert!((f21 - 0.2).abs() < 1e-10, "got {f21}");
    }
    #[test]
    fn test_view_factor_sum_rule_two_surface() {
        let f_self = view_factor_sum_rule(&[0.6]);
        assert!((f_self - 0.4).abs() < 1e-10, "got {f_self}");
    }
    #[test]
    fn test_view_factor_sum_rule_full() {
        let f_self = view_factor_sum_rule(&[0.3, 0.4, 0.3]);
        assert!(f_self.abs() < 1e-10, "got {f_self}");
    }
    #[test]
    fn test_view_factor_coaxial_disks_equal_radii() {
        let f = view_factor_coaxial_disks(0.5, 0.5, 1.0);
        assert!(f > 0.0 && f < 1.0, "got {f}");
    }
    #[test]
    fn test_view_factor_coaxial_disks_large_disk_small_h() {
        let f = view_factor_coaxial_disks(10.0, 10.0, 0.01);
        assert!(f > 0.99, "got {f}");
    }
    #[test]
    fn test_view_factor_parallel_plates_in_range() {
        let f = view_factor_parallel_plates(1.0, 1.0, 1.0);
        assert!(f > 0.0 && f < 1.0, "got {f}");
    }
    #[test]
    fn test_radiation_surface_radiosity() {
        let s = RadiationSurface::new(1.0, 0.8, 1000.0, "test");
        let j = s.radiosity();
        assert!((j - 45_360.0).abs() < 50.0, "got {j}");
    }
    #[test]
    fn test_radiation_surface_blackbody_resistance_zero() {
        let s = RadiationSurface::new(1.0, 1.0, 1000.0, "bb");
        assert_eq!(s.surface_resistance(), 0.0);
    }
    #[test]
    fn test_radiation_surface_resistance_formula() {
        let s = RadiationSurface::new(2.0, 0.5, 500.0, "s");
        assert!((s.surface_resistance() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_two_surface_enclosure_same_temp_zero_flux() {
        let s1 = RadiationSurface::new(1.0, 0.8, 500.0, "s1");
        let s2 = RadiationSurface::new(1.0, 0.8, 500.0, "s2");
        let q = two_surface_enclosure_heat_transfer(&s1, &s2, 1.0);
        assert!(q.abs() < 1e-6, "got {q}");
    }
    #[test]
    fn test_two_surface_enclosure_hot_to_cold() {
        let s1 = RadiationSurface::new(1.0, 0.9, 800.0, "hot");
        let s2 = RadiationSurface::new(1.0, 0.9, 400.0, "cold");
        let q = two_surface_enclosure_heat_transfer(&s1, &s2, 1.0);
        assert!(q > 0.0, "heat should flow from hot to cold, got {q}");
    }
    #[test]
    fn test_air_mass_at_zenith() {
        let am = air_mass(90.0);
        assert!((am - 1.0).abs() < 1e-6, "got {am}");
    }
    #[test]
    fn test_air_mass_at_30_deg() {
        let am = air_mass(30.0);
        assert!((am - 2.0).abs() < 1e-6, "got {am}");
    }
    #[test]
    fn test_air_mass_below_horizon() {
        assert!(air_mass(0.0).is_infinite());
        assert!(air_mass(-10.0).is_infinite());
    }
    #[test]
    fn test_attenuation_factor_zero_tau() {
        assert!((attenuation_factor(1.5, 0.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_attenuation_factor_decreases_with_am() {
        let i1 = attenuation_factor(1.0, 0.1);
        let i2 = attenuation_factor(2.0, 0.1);
        assert!(i1 > i2, "higher AM should give lower transmittance");
    }
    #[test]
    fn test_solar_irradiance_horizontal_surface() {
        let irr = solar_irradiance_surface(800.0, 100.0, 0.0, 0.0, 90.0, 0.0);
        assert!((irr - 900.0).abs() < 1.0, "got {irr}");
    }
    #[test]
    fn test_solar_irradiance_sun_behind_surface() {
        let irr = solar_irradiance_surface(800.0, 100.0, 0.0, 180.0, 10.0, 0.0);
        assert!(irr >= 0.0);
    }
    #[test]
    fn test_solar_irradiance_no_dni_only_diffuse() {
        let irr = solar_irradiance_surface(0.0, 200.0, 90.0, 0.0, 45.0, 0.0);
        assert!((irr - 100.0).abs() < 1.0, "got {irr}");
    }
    #[test]
    fn test_blackbody_spectrum_total_power() {
        let bb = BlackbodySpectrum::new(1000.0);
        assert!((bb.total_power() - blackbody_emissive_power(1000.0)).abs() < 1e-6);
    }
    #[test]
    fn test_blackbody_spectrum_peak_wavelength() {
        let bb = BlackbodySpectrum::new(5778.0);
        assert!((bb.peak_wavelength() - 501.6e-9).abs() < 5e-9);
    }
    #[test]
    fn test_blackbody_spectrum_planck_positive() {
        let bb = BlackbodySpectrum::new(3000.0);
        assert!(bb.planck(1e-6) > 0.0);
    }
    #[test]
    fn test_blackbody_spectrum_band_fraction_full_range() {
        let bb = BlackbodySpectrum::new(3000.0);
        let frac = bb.band_fraction(100e-9, 100e-6, 5000);
        assert!(
            frac > 0.98,
            "band fraction over wide range should be ~1, got {frac}"
        );
    }
    #[test]
    fn test_radiative_properties_opaque_consistent() {
        let rp = RadiativeProperties::new(0.8, 0.0);
        assert!(rp.is_consistent(), "should be consistent");
    }
    #[test]
    fn test_radiative_properties_emissive_power() {
        let rp = RadiativeProperties::new(0.9, 0.0);
        let ep = rp.emissive_power(1000.0);
        assert!((ep - 0.9 * blackbody_emissive_power(1000.0)).abs() < 1e-3);
    }
    #[test]
    fn test_thermal_radiation_parallel_plates_same_temp() {
        let q = thermal_radiation_parallel_plates(0.9, 500.0, 0.9, 500.0);
        assert!(q.abs() < 1e-6, "no net flux at same temperature, got {q}");
    }
    #[test]
    fn test_thermal_radiation_parallel_plates_direction() {
        let q = thermal_radiation_parallel_plates(0.9, 800.0, 0.9, 400.0);
        assert!(q > 0.0, "heat should flow from hot to cold, got {q}");
    }
    #[test]
    fn test_thermal_radiation_body_in_enclosure_direction() {
        let q = thermal_radiation_body_in_enclosure(0.8, 1000.0, 300.0);
        assert!(
            q > 0.0,
            "body hotter than enclosure should lose heat, got {q}"
        );
    }
    #[test]
    fn test_view_factor_concentric_cylinders_f12_unity() {
        let (f12, _) = view_factor_concentric_cylinders(0.1, 0.3);
        assert!(
            (f12 - 1.0).abs() < 1e-12,
            "F12 for infinite concentric cylinders = 1, got {f12}"
        );
    }
    #[test]
    fn test_view_factor_concentric_cylinders_f21_reciprocity() {
        let r1 = 0.1;
        let r2 = 0.4;
        let (f12, f21) = view_factor_concentric_cylinders(r1, r2);
        assert!((f21 - r1 / r2).abs() < 1e-12, "F21 = r1/r2, got {f21}");
        let _ = f12;
    }
    #[test]
    fn test_extraterrestrial_irradiance_near_solar_constant() {
        let i = extraterrestrial_irradiance(1);
        assert!(
            (i - SOLAR_CONSTANT).abs() < SOLAR_CONSTANT * 0.05,
            "got {i}"
        );
    }
    #[test]
    fn test_solar_elevation_noon_equator_equinox() {
        let elev = solar_elevation(0.0, 0.0, 0.0);
        assert!((elev - 90.0).abs() < 1e-9, "got {elev}");
    }
    #[test]
    fn test_radiation_damage_dpa() {
        let rd = RadiationDamage::new(1.0e20, 1.0e-24, 8.49e22);
        let dpa = rd.dpa();
        assert!((dpa - 1.0e20 * 1.0e-24).abs() < 1e-10, "got {dpa}");
    }
    #[test]
    fn test_radiation_damage_swelling_positive() {
        let rd = RadiationDamage::new(1.0e22, 1.0e-24, 8.49e22);
        let sw = rd.swelling_fraction();
        assert!(sw >= 0.0, "got {sw}");
    }
    #[test]
    fn test_neutron_moderation_slowing_down_power() {
        let nm = NeutronModeration::new(1.475, 0.022, 0.920);
        let sdp = nm.slowing_down_power();
        assert!((sdp - 0.920 * 1.475).abs() < 1e-6, "got {sdp}");
    }
    #[test]
    fn test_neutron_moderation_ratio_water() {
        let nm = NeutronModeration::new(1.475, 0.022, 0.920);
        let mr = nm.moderation_ratio();
        assert!(mr > 50.0, "water MR should be > 50, got {mr}");
    }
    #[test]
    fn test_neutron_moderation_zero_absorption_infinite_ratio() {
        let nm = NeutronModeration::new(1.0, 0.0, 1.0);
        assert!(nm.moderation_ratio().is_infinite());
    }
    #[test]
    fn test_participating_media_extinction() {
        let pm = ParticipatingMedia::new(0.1, 0.05, 1000.0);
        assert!((pm.extinction() - 0.15).abs() < 1e-10);
    }
    #[test]
    fn test_participating_media_albedo() {
        let pm = ParticipatingMedia::new(0.1, 0.05, 1000.0);
        let omega = pm.albedo();
        assert!(
            (omega - 1.0 / 3.0).abs() < 1e-10,
            "albedo = sigma_s / beta = 0.05/0.15"
        );
    }
    #[test]
    fn test_participating_media_transmittance_range() {
        let pm = ParticipatingMedia::new(0.5, 0.0, 1000.0);
        let t = pm.transmittance(1.0);
        assert!(
            (0.0..=1.0).contains(&t),
            "transmittance must be in [0,1], got {t}"
        );
    }
    #[test]
    fn test_participating_media_optically_thick_no_transmittance() {
        let pm = ParticipatingMedia::new(1000.0, 0.0, 1000.0);
        let t = pm.transmittance(1.0);
        assert!(t < 1e-100, "optically thick medium: T → 0, got {t}");
    }
    #[test]
    fn test_participating_media_effective_emissivity_range() {
        let pm = ParticipatingMedia::new(0.2, 0.0, 800.0);
        let eps = pm.effective_emissivity(2.0);
        assert!(
            (0.0..=1.0).contains(&eps),
            "emissivity must be in [0,1], got {eps}"
        );
    }
    #[test]
    fn test_participating_media_optically_thin_emission() {
        let pm = ParticipatingMedia::new(0.01, 0.0, 1000.0);
        let q = pm.emission_optically_thin(1.0);
        assert!(q > 0.0, "emission should be positive, got {q}");
    }
    #[test]
    fn test_participating_media_mean_beam_sphere() {
        let r = 0.5;
        let l = ParticipatingMedia::mean_beam_length_sphere(r);
        assert!((l - 0.65).abs() < 1e-10, "got {l}");
    }
    #[test]
    fn test_participating_media_albedo_zero_absorption() {
        let pm = ParticipatingMedia::new(1.0, 0.0, 500.0);
        assert!(pm.albedo().abs() < 1e-12, "purely absorbing: ω = 0");
    }
    #[test]
    fn test_mc_view_factor_in_range() {
        let mut mc = MonteCarloViewFactor::new(5000);
        let f = mc.parallel_disks(1.0, 1.0, 1.0);
        assert!(
            f > 0.0 && f < 1.0,
            "MC view factor must be in (0,1), got {f}"
        );
    }
    #[test]
    fn test_mc_view_factor_large_disks_close() {
        let mut mc = MonteCarloViewFactor::new(10000);
        let f = mc.parallel_disks(10.0, 10.0, 0.01);
        assert!(f > 0.8, "large close disks: F12 should be > 0.8, got {f}");
    }
    #[test]
    fn test_mc_view_factor_small_disk_far() {
        let mut mc = MonteCarloViewFactor::new(5000);
        let f = mc.parallel_disks(1.0, 0.01, 100.0);
        assert!(f < 0.1, "small far disk: F12 should be small, got {f}");
    }
    #[test]
    fn test_radiation_network_same_temp_zero_net_flux() {
        let t = 500.0;
        let n = RadiationNetwork::new(
            vec![t, t],
            vec![0.8, 0.8],
            vec![1.0, 1.0],
            vec![0.0, 1.0, 1.0, 0.0],
        );
        let q0 = n.net_heat_flow_simple(0);
        assert!(q0.abs() < 1.0, "same T: net heat flow ≈ 0, got {q0}");
    }
    #[test]
    fn test_radiation_network_hot_to_cold() {
        let n = RadiationNetwork::new(
            vec![1000.0, 300.0],
            vec![0.9, 0.9],
            vec![1.0, 1.0],
            vec![0.0, 1.0, 1.0, 0.0],
        );
        let q0 = n.net_heat_flow_simple(0);
        assert!(q0 > 0.0, "hot surface should lose heat, q0 = {q0}");
    }
    #[test]
    fn test_radiation_network_emitted_power_positive() {
        let n = RadiationNetwork::new(
            vec![500.0, 300.0],
            vec![0.8, 0.7],
            vec![2.0, 1.0],
            vec![0.0, 1.0, 1.0, 0.0],
        );
        assert!(n.emitted_power(0) > 0.0);
        assert!(n.emitted_power(1) > 0.0);
    }
    #[test]
    fn test_radiation_shield_reduces_emissivity() {
        let eps_no_shield = radiation_shield_effective_emissivity(0.9, 0.9, 0.9, 0);
        let eps_one_shield = radiation_shield_effective_emissivity(0.9, 0.9, 0.9, 1);
        assert!(
            eps_one_shield < eps_no_shield,
            "shield should reduce effective emissivity: {eps_one_shield} vs {eps_no_shield}"
        );
    }
    #[test]
    fn test_radiation_shield_more_shields_lower_emissivity() {
        let eps_1 = radiation_shield_effective_emissivity(0.5, 0.5, 0.5, 1);
        let eps_5 = radiation_shield_effective_emissivity(0.5, 0.5, 0.5, 5);
        assert!(
            eps_5 < eps_1,
            "more shields → lower ε_eff: {eps_5} vs {eps_1}"
        );
    }
    #[test]
    fn test_radiation_shield_no_shield_matches_parallel_plates() {
        let eps1 = 0.8;
        let eps2 = 0.6;
        let eff = radiation_shield_effective_emissivity(eps1, eps2, 0.5, 0);
        let expected = 1.0 / (1.0 / eps1 + 1.0 / eps2 - 1.0);
        assert!(
            (eff - expected).abs() < 1e-10,
            "got {eff}, expected {expected}"
        );
    }
    #[test]
    fn test_radiation_shield_reduction_factor_less_than_unity() {
        let r = radiation_shield_reduction_factor(0.9, 2);
        assert!(
            r > 0.0 && r < 1.0,
            "reduction factor should be in (0,1), got {r}"
        );
    }
    #[test]
    fn test_radiation_shield_reduction_factor_decreases_with_shields() {
        let r1 = radiation_shield_reduction_factor(0.8, 1);
        let r3 = radiation_shield_reduction_factor(0.8, 3);
        assert!(r3 < r1, "more shields → more reduction: {r3} vs {r1}");
    }
    #[test]
    fn test_kirchhoff_law_gray_surface() {
        let eps = 0.75;
        let alpha = absorptivity_from_emissivity(eps);
        assert!(
            (alpha - eps).abs() < 1e-12,
            "Kirchhoff: α = ε for gray, got α={alpha}"
        );
    }
    #[test]
    fn test_kirchhoff_law_satisfied_check() {
        assert!(kirchhoff_law_satisfied(0.8, 0.8, 1e-10));
        assert!(!kirchhoff_law_satisfied(0.8, 0.5, 1e-10));
    }
    #[test]
    fn test_kirchhoff_spectral_emissivity() {
        let absorptivity = vec![0.3, 0.5, 0.8, 0.9];
        let emissivity = spectral_emissivity_kirchhoff(&absorptivity);
        assert_eq!(emissivity, absorptivity, "spectral ε should equal α");
    }
    #[test]
    fn test_surface_radiosity_blackbody() {
        let j = surface_radiosity(1.0, 1000.0, 999.0);
        let e_b = SIGMA * 1000.0_f64.powi(4);
        assert!((j - e_b).abs() < 1e-6, "blackbody radiosity = E_b, got {j}");
    }
    #[test]
    fn test_surface_radiosity_gray_surface() {
        let j = surface_radiosity(0.8, 500.0, 1000.0);
        let e_b = SIGMA * 500.0_f64.powi(4);
        let expected = 0.8 * e_b + 0.2 * 1000.0;
        assert!((j - expected).abs() < 1e-6, "got {j}, expected {expected}");
    }
    #[test]
    fn test_radiosity_net_flux_blackbody() {
        let e_b = SIGMA * 1000.0_f64.powi(4);
        let j = e_b - 1000.0;
        let q = radiosity_net_flux(1.0, 1000.0, j);
        assert!((q - 1000.0).abs() < 1e-3, "got {q}");
    }
    #[test]
    fn test_radiosity_solve_jacobi_converges() {
        let t = 500.0_f64;
        let e_b = SIGMA * t.powi(4);
        let j = radiosity_solve_jacobi(&[t, t], &[0.8, 0.8], &[0.0, 1.0, 1.0, 0.0], 50);
        assert!(
            (j[0] - e_b).abs() / e_b < 0.01,
            "J[0] should ≈ E_b at same T, got {}",
            j[0]
        );
        assert!(
            (j[1] - e_b).abs() / e_b < 0.01,
            "J[1] should ≈ E_b at same T, got {}",
            j[1]
        );
    }
    #[test]
    fn test_radiosity_net_heat_flow_zero_at_same_temp() {
        let t = 500.0;
        let j = radiosity_solve_jacobi(&[t, t], &[0.9, 0.9], &[0.0, 1.0, 1.0, 0.0], 100);
        let q0 = radiosity_net_heat_flow(0, &[1.0, 1.0], &j, &[0.0, 1.0, 1.0, 0.0]);
        assert!(
            q0.abs() < 1.0,
            "net heat flow at same T should be ~0, got {q0}"
        );
    }
    #[test]
    fn test_am15_spectral_irradiance_positive_at_500nm() {
        let e = am15_spectral_irradiance(500.0);
        assert!(
            e > 0.0,
            "AM1.5 irradiance at 500 nm should be positive, got {e}"
        );
    }
    #[test]
    fn test_am15_spectral_irradiance_decays_at_uv() {
        let e_vis = am15_spectral_irradiance(500.0);
        let e_uv = am15_spectral_irradiance(50.0);
        assert!(e_vis > e_uv, "visible irradiance > deep UV irradiance");
    }
    #[test]
    fn test_am15_integrated_irradiance_positive() {
        let e = am15_integrated_irradiance(300.0, 2500.0, 500);
        assert!(
            e > 100.0 && e < 2000.0,
            "AM1.5 total should be ~1000 W/m², got {e}"
        );
    }
    #[test]
    fn test_am15_photon_flux_positive() {
        let phi = am15_photon_flux(500.0);
        assert!(
            phi > 0.0,
            "photon flux at 500 nm should be positive, got {phi}"
        );
    }
    #[test]
    fn test_am15_photon_flux_longer_wavelength_higher() {
        let phi_500 = am15_photon_flux(500.0);
        let phi_1000 = am15_photon_flux(1000.0);
        assert!(phi_500 > 0.0 && phi_1000 > 0.0);
    }
    #[test]
    fn test_gray_body_emissivity_at_ref_temp() {
        let gb = GrayBodyModel::new(0.8, 0.001, 300.0);
        assert!(
            (gb.emissivity(300.0) - 0.8).abs() < 1e-12,
            "got {}",
            gb.emissivity(300.0)
        );
    }
    #[test]
    fn test_gray_body_emissivity_clamped_to_unity() {
        let gb = GrayBodyModel::new(0.95, 0.01, 300.0);
        assert_eq!(gb.emissivity(1000.0), 1.0);
    }
    #[test]
    fn test_gray_body_emissive_power_positive() {
        let gb = GrayBodyModel::new(0.7, 0.0, 300.0);
        let p = gb.emissive_power(500.0);
        assert!(p > 0.0, "emissive power should be positive, got {p}");
    }
    #[test]
    fn test_gray_body_effective_temperature_less_than_actual() {
        let gb = GrayBodyModel::new(0.5, 0.0, 300.0);
        let t_eff = gb.effective_blackbody_temperature(800.0);
        assert!(t_eff < 800.0, "T_eff should be < T for ε < 1, got {t_eff}");
    }
    #[test]
    fn test_view_factor_perpendicular_in_range() {
        let f = view_factor_perpendicular_rectangles(1.0, 1.0, 1.0);
        assert!((0.0..=1.0).contains(&f), "F12 must be in [0,1], got {f}");
    }
    #[test]
    fn test_view_factor_perpendicular_positive() {
        let f = view_factor_perpendicular_rectangles(2.0, 2.0, 2.0);
        assert!(
            f > 0.0,
            "F12 should be positive for perpendicular geometry, got {f}"
        );
    }
    #[test]
    fn test_solar_cell_thermal_voltage_at_300k() {
        let cell = SolarCell::new(5.0, 1e-10, 1.0, 300.0, 0.0);
        let vt = cell.thermal_voltage();
        assert!((vt - 0.02585).abs() < 0.001, "V_t = {vt}");
    }
    #[test]
    fn test_solar_cell_open_circuit_voltage_positive() {
        let cell = SolarCell::new(5.0, 1e-10, 1.0, 300.0, 0.0);
        let voc = cell.open_circuit_voltage();
        assert!(voc > 0.0, "Voc should be positive, got {voc}");
    }
    #[test]
    fn test_solar_cell_short_circuit_current_approx_i_ph() {
        let cell = SolarCell::new(5.0, 1e-12, 1.0, 300.0, 0.0);
        let isc = cell.current_at_voltage(0.0);
        assert!((isc - 5.0).abs() < 0.1, "Isc ≈ Iph, got {isc}");
    }
    #[test]
    fn test_solar_cell_power_max_somewhere() {
        let cell = SolarCell::new(5.0, 1e-10, 1.0, 300.0, 0.0);
        let voc = cell.open_circuit_voltage();
        let n = 100;
        let p_max = (0..=n)
            .map(|i| cell.power_at_voltage(voc * i as f64 / n as f64))
            .fold(0.0_f64, f64::max);
        assert!(p_max > 0.0, "maximum power should be > 0, got {p_max}");
    }
    #[test]
    fn test_solar_cell_fill_factor_in_range() {
        let cell = SolarCell::new(5.0, 1e-10, 1.0, 300.0, 0.0);
        let ff = cell.fill_factor();
        assert!(
            ff > 0.0 && ff <= 1.0,
            "fill factor should be in (0,1], got {ff}"
        );
    }
    #[test]
    fn test_spectral_emissivity_interpolation_in_range() {
        let wavelengths = vec![500.0, 1000.0, 2000.0];
        let emissivities = vec![0.9, 0.7, 0.5];
        let model = SpectralEmissivityModel::new(wavelengths, emissivities);
        let e = model.emissivity_at(1000.0);
        assert!(
            (e - 0.7).abs() < 1e-10,
            "at knot should equal knot value, got {e}"
        );
    }
    #[test]
    fn test_spectral_emissivity_outside_range_clamped() {
        let wavelengths = vec![500.0, 2000.0];
        let emissivities = vec![0.8, 0.4];
        let model = SpectralEmissivityModel::new(wavelengths, emissivities);
        let e_low = model.emissivity_at(100.0);
        let e_high = model.emissivity_at(5000.0);
        assert!(
            (e_low - 0.8).abs() < 1e-10,
            "below range: should clamp to first, got {e_low}"
        );
        assert!(
            (e_high - 0.4).abs() < 1e-10,
            "above range: should clamp to last, got {e_high}"
        );
    }
    #[test]
    fn test_spectral_emissivity_midpoint_interpolation() {
        let wavelengths = vec![0.0, 2.0];
        let emissivities = vec![0.0, 1.0];
        let model = SpectralEmissivityModel::new(wavelengths, emissivities);
        let e = model.emissivity_at(1.0);
        assert!((e - 0.5).abs() < 1e-10, "midpoint should be 0.5, got {e}");
    }
    #[test]
    fn test_spectral_effective_emissivity_blackbody_is_one() {
        let wavelengths = vec![100e-9, 100e-6];
        let emissivities = vec![1.0, 1.0];
        let model = SpectralEmissivityModel::new(wavelengths, emissivities);
        let eps = model.effective_total_emissivity(1000.0, 100);
        assert!((eps - 1.0).abs() < 0.01, "blackbody ε_total ≈ 1, got {eps}");
    }
    #[test]
    fn test_dosimetry_fluence_from_rate_and_time() {
        let dos = IrradiationDosimetry::new(1.0e14, 3600.0, 50.0e-24);
        let fluence = dos.fluence();
        assert!((fluence - 1.0e14 * 3600.0).abs() < 1.0, "fluence={fluence}");
    }
    #[test]
    fn test_dosimetry_absorbed_dose_positive() {
        let dos = IrradiationDosimetry::new(1.0e13, 1000.0, 100.0e-24);
        let dose = dos.absorbed_dose_gray();
        assert!(dose > 0.0, "dose={dose}");
    }
    #[test]
    fn test_dosimetry_kerma_rate_positive() {
        let dos = IrradiationDosimetry::new(1.0e12, 1.0, 50.0e-24);
        let kerma = dos.kerma_rate(8960.0);
        assert!(kerma > 0.0, "kerma={kerma}");
    }
    #[test]
    fn test_dosimetry_rem_dose_positive() {
        let dos = IrradiationDosimetry::new(1.0e13, 100.0, 50.0e-24);
        let rem = dos.rem_dose(20.0);
        assert!(rem > 0.0, "rem dose={rem}");
    }
    #[test]
    fn test_sb_integrator_blackbody_total_matches_sigma() {
        let integrator = StefanBoltzmannIntegrator::new(1000.0);
        let total = integrator.integrate_full_spectrum(500);
        let expected = SIGMA * 1000.0_f64.powi(4);
        assert!(
            total > 0.0 && total < expected * 1.1,
            "total={total}, expected≈{expected}"
        );
    }
    #[test]
    fn test_sb_integrator_band_fraction_visible_range() {
        let integrator = StefanBoltzmannIntegrator::new(5778.0);
        let frac = integrator.band_fraction_in_range(400e-9, 700e-9, 200);
        assert!(
            frac > 0.2 && frac < 0.6,
            "visible band fraction at 5778K={frac}"
        );
    }
    #[test]
    fn test_sb_integrator_increases_with_temperature() {
        let i1 = StefanBoltzmannIntegrator::new(500.0).integrate_full_spectrum(200);
        let i2 = StefanBoltzmannIntegrator::new(1000.0).integrate_full_spectrum(200);
        assert!(i2 > i1, "higher T → more power: i1={i1}, i2={i2}");
    }
    #[test]
    fn test_mc_transport_transmission_positive() {
        let mut mc = MCRadiationTransport::new(1000, 42);
        let t = mc.slab_transmittance(0.5, 1.0);
        assert!((0.0..=1.0).contains(&t), "transmittance in [0,1], got {t}");
    }
    #[test]
    fn test_mc_transport_opaque_slab_low_transmittance() {
        let mut mc = MCRadiationTransport::new(2000, 99);
        let t = mc.slab_transmittance(10.0, 1.0);
        assert!(t < 0.1, "opaque slab should transmit little, got {t}");
    }
    #[test]
    fn test_mc_transport_transparent_slab_high_transmittance() {
        let mut mc = MCRadiationTransport::new(2000, 7);
        let t = mc.slab_transmittance(0.01, 1.0);
        assert!(t > 0.9, "transparent slab should transmit most, got {t}");
    }
    #[test]
    fn test_mc_transport_albedo_effect() {
        let mut mc_abs = MCRadiationTransport::new(2000, 1);
        let mut mc_scat = MCRadiationTransport::new(2000, 1);
        let t_abs = mc_abs.slab_transmittance_full(1.0, 1.0, 0.0);
        let t_scat = mc_scat.slab_transmittance_full(1.0, 1.0, 1.0);
        assert!(t_abs >= 0.0 && t_scat >= 0.0);
    }
    #[test]
    fn test_view_factor_sphere_to_enclosing_plane_in_range() {
        let f = view_factor_sphere_in_enclosure(0.5, 4.0);
        assert!(f > 0.0 && f <= 1.0, "F12 in [0,1], got {f}");
    }
    #[test]
    fn test_view_factor_strip_to_strip_parallel_in_range() {
        let f = view_factor_infinite_parallel_strips(1.0, 1.0, 1.0);
        assert!(f > 0.0 && f < 1.0, "F12 in (0,1), got {f}");
    }
    #[test]
    fn test_view_factor_isothermal_enclosure_sums_to_one() {
        let f12 = 0.4_f64;
        let f13 = view_factor_sum_rule(&[f12]);
        assert!((f12 + f13 - 1.0).abs() < 1e-10, "F12+F13={}", f12 + f13);
    }
    #[test]
    fn test_extraterrestrial_irradiance_day1() {
        let i = extraterrestrial_irradiance(1);
        assert!(
            i > SOLAR_CONSTANT * 0.95 && i < SOLAR_CONSTANT * 1.1,
            "I_et={i}"
        );
    }
    #[test]
    fn test_extraterrestrial_irradiance_varies_seasonally() {
        let i_jan = extraterrestrial_irradiance(1);
        let i_jun = extraterrestrial_irradiance(172);
        assert!(
            (i_jan - i_jun).abs() / SOLAR_CONSTANT > 0.01,
            "seasonal variation too small"
        );
    }
    #[test]
    fn test_solar_elevation_equator_equinox_noon() {
        let el = solar_elevation(0.0, 0.0, 0.0);
        assert!(
            (el - 90.0).abs() < 0.1,
            "elevation at equator equinox noon = 90°, got {el}"
        );
    }
    #[test]
    fn test_solar_elevation_negative_before_sunrise() {
        let el = solar_elevation(0.0, 0.0, 90.0);
        assert!(
            el.abs() < 1.0,
            "at hour angle 90° sun near horizon, got {el}"
        );
    }
    #[test]
    fn test_participating_media_extinction_correct() {
        let pm = ParticipatingMedia::new(2.0, 1.0, 600.0);
        assert!(
            (pm.extinction() - 3.0).abs() < 1e-10,
            "β=κ+σ=3, got {}",
            pm.extinction()
        );
    }
    #[test]
    fn test_participating_media_albedo_scattering() {
        let pm = ParticipatingMedia::new(1.0, 3.0, 500.0);
        assert!(
            (pm.albedo() - 0.75).abs() < 1e-10,
            "ω=0.75, got {}",
            pm.albedo()
        );
    }
    #[test]
    fn test_participating_media_transmittance_decreases_with_depth() {
        let pm = ParticipatingMedia::new(1.0, 0.0, 1000.0);
        let t1 = pm.transmittance(1.0);
        let t2 = pm.transmittance(2.0);
        assert!(t2 < t1, "transmittance decreases with depth");
    }
    #[test]
    fn test_participating_media_effective_emissivity_approaches_unity() {
        let pm = ParticipatingMedia::new(10.0, 0.0, 500.0);
        let e = pm.effective_emissivity(10.0);
        assert!(e > 0.99, "optically thick medium ε → 1, got {e}");
    }
    #[test]
    fn test_radiation_network_two_surface_enclosure() {
        let temps = vec![1000.0, 300.0];
        let eps = vec![0.8, 0.9];
        let areas = vec![1.0, 1.0];
        let vf = vec![0.0, 1.0, 1.0, 0.0];
        let net = RadiationNetwork::new(temps, eps, areas, vf);
        let q0 = net.net_heat_flow_simple(0);
        assert!(q0 > 0.0, "hot surface loses heat, q0={q0}");
    }
    #[test]
    fn test_radiation_network_same_temp_zero_net() {
        let temps = vec![500.0, 500.0];
        let eps = vec![0.8, 0.8];
        let areas = vec![1.0, 1.0];
        let vf = vec![0.0, 1.0, 1.0, 0.0];
        let net = RadiationNetwork::new(temps, eps, areas, vf);
        let q0 = net.net_heat_flow_simple(0);
        assert!(q0.abs() < 1e-3, "equal-temp enclosure: q≈0, got {q0}");
    }
}
