//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{GtnParams, PhaseFieldParams};

/// Compute the crack driving force H⁺ = max(ψ₀, H_prev).
///
/// In the AT2 model the history variable H stores the maximum tensile strain
/// energy density ever reached, ensuring irreversibility of crack growth.
#[allow(dead_code)]
pub fn crack_driving_force(strain_energy_density: f64, h_prev: f64) -> f64 {
    strain_energy_density.max(h_prev)
}
/// Degradation function g(d) = (1-d)² + k_res, where k_res = 1e-6.
///
/// Multiplies the elastic stiffness to model loss of load-carrying capacity.
/// The small residual k_res prevents numerical singularity at d=1.
#[allow(dead_code)]
pub fn degradation_function(d: f64) -> f64 {
    pub(super) const K_RES: f64 = 1.0e-6;
    (1.0 - d).powi(2) + K_RES
}
/// Phase field residual for the AT2 model.
///
/// R(d) = (Gc/l0) * d − 2*(1−d)*H
///
/// Solving R=0 gives the equilibrium phase field value at a given history H.
#[allow(dead_code)]
pub fn phase_field_residual(d: f64, h: f64, params: &PhaseFieldParams) -> f64 {
    (params.gc / params.l0) * d - 2.0 * (1.0 - d) * h
}
/// Update the phase field d by solving the AT2 model in closed form.
///
/// Setting the residual R(d) = 0 and solving:
///   d = 2*H / (Gc/l0 + 2*H)
///
/// The result is clamped to \[0, 1\] and irreversibility is enforced (d only
/// increases).
#[allow(dead_code)]
pub fn update_phase_field(d_prev: f64, h: f64, params: &PhaseFieldParams) -> f64 {
    let alpha = params.gc / params.l0;
    let d_new = 2.0 * h / (alpha + 2.0 * h);
    d_new.clamp(d_prev, 1.0)
}
/// Effective (degraded) Young's modulus: E_eff = g(d) * E0.
#[allow(dead_code)]
pub fn effective_modulus(d: f64, e0: f64) -> f64 {
    degradation_function(d) * e0
}
/// Approximate crack surface energy from a 1-D phase field array.
///
/// Γ ≈ Σᵢ \[ dᵢ²/(2*l0) + (l0/2)*(∂d/∂x)² \] * dx
///
/// Uses central differences for interior nodes and one-sided differences at
/// boundaries.
#[allow(dead_code)]
pub fn crack_length_from_field(d_field: &[f64], dx: f64, l0: f64) -> f64 {
    let n = d_field.len();
    if n == 0 {
        return 0.0;
    }
    let mut gamma = 0.0;
    for i in 0..n {
        let grad = if n == 1 {
            0.0
        } else if i == 0 {
            (d_field[1] - d_field[0]) / dx
        } else if i == n - 1 {
            (d_field[n - 1] - d_field[n - 2]) / dx
        } else {
            (d_field[i + 1] - d_field[i - 1]) / (2.0 * dx)
        };
        gamma += (d_field[i].powi(2) / (2.0 * l0) + (l0 / 2.0) * grad.powi(2)) * dx;
    }
    gamma
}
/// Mode I SIF for a through-thickness center crack in an infinite plate.
///
/// K_I = σ * √(π * a)
///
/// where σ is the applied stress and a is the half-crack length.
#[allow(dead_code)]
pub fn sif_center_crack_infinite(stress: f64, half_crack_length: f64) -> f64 {
    stress * (PI * half_crack_length).sqrt()
}
/// Mode I SIF for a through-thickness center crack in a finite-width plate.
///
/// K_I = σ * √(π * a) * F(a/W)
///
/// where W is the half-width and F(a/W) is the finite-width correction factor:
/// F(a/W) ≈ (1 - 0.025*(a/W)^2 + 0.06*(a/W)^4) / cos(π*a/(2W))^(1/2)
/// (Feddersen/Isida approximation)
#[allow(dead_code)]
pub fn sif_center_crack_finite_width(stress: f64, half_crack: f64, half_width: f64) -> f64 {
    let ratio = half_crack / half_width;
    let correction =
        (1.0 - 0.025 * ratio.powi(2) + 0.06 * ratio.powi(4)) / (PI * ratio / 2.0).cos().sqrt();
    stress * (PI * half_crack).sqrt() * correction
}
/// Mode I SIF for a single edge crack in a semi-infinite plate.
///
/// K_I = 1.12 * σ * √(π * a)
#[allow(dead_code)]
pub fn sif_single_edge_crack(stress: f64, crack_length: f64) -> f64 {
    1.12 * stress * (PI * crack_length).sqrt()
}
/// Mode I SIF for a single edge crack in a finite-width plate.
///
/// K_I = σ * √(π * a) * F(a/W)
///
/// F(a/W) = 1.12 - 0.231*(a/W) + 10.55*(a/W)^2 - 21.71*(a/W)^3 + 30.38*(a/W)^4
/// (Brown and Srawley polynomial fit for a/W < 0.6)
#[allow(dead_code)]
pub fn sif_edge_crack_finite_width(stress: f64, crack_length: f64, width: f64) -> f64 {
    let r = crack_length / width;
    let f = 1.12 - 0.231 * r + 10.55 * r.powi(2) - 21.71 * r.powi(3) + 30.38 * r.powi(4);
    stress * (PI * crack_length).sqrt() * f
}
/// Mode I SIF for a double edge crack in a finite-width plate.
///
/// K_I = σ * √(π * a) * F(a/W)
///
/// F(a/W) = 1.12 + 0.203*(a/W) - 1.197*(a/W)^2 + 1.93*(a/W)^3
/// (Tada, Paris and Irwin approximation)
#[allow(dead_code)]
pub fn sif_double_edge_crack(stress: f64, crack_length: f64, half_width: f64) -> f64 {
    let r = crack_length / half_width;
    let f = 1.12 + 0.203 * r - 1.197 * r.powi(2) + 1.93 * r.powi(3);
    stress * (PI * crack_length).sqrt() * f
}
/// Mode I SIF for a penny-shaped (circular) crack in an infinite body.
///
/// K_I = (2/π) * σ * √(π * a)
///
/// where a is the crack radius.
#[allow(dead_code)]
pub fn sif_penny_shaped_crack(stress: f64, crack_radius: f64) -> f64 {
    (2.0 / PI) * stress * (PI * crack_radius).sqrt()
}
/// Mode I SIF for a compact tension (CT) specimen (ASTM E399).
///
/// K_I = (P / (B * √W)) * f(a/W)
///
/// f(a/W) = (2 + α) * (0.886 + 4.64*α - 13.32*α^2 + 14.72*α^3 - 5.6*α^4) / (1-α)^(3/2)
/// where α = a/W.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn sif_compact_tension(load: f64, thickness: f64, width: f64, crack_length: f64) -> f64 {
    let alpha = crack_length / width;
    let f = (2.0 + alpha)
        * (0.886 + 4.64 * alpha - 13.32 * alpha.powi(2) + 14.72 * alpha.powi(3)
            - 5.6 * alpha.powi(4))
        / (1.0 - alpha).powf(1.5);
    (load / (thickness * width.sqrt())) * f
}
/// Evaluate the Bueckner weight function for an edge crack.
///
/// m(x, a) = (2 / √(2π(a-x))) * (1 + M1*(1-x/a)^(1/2) + M2*(1-x/a) + M3*(1-x/a)^(3/2))
///
/// where M1, M2, M3 are geometry-dependent coefficients.
/// For a standard edge crack: M1 ≈ 0.6147, M2 ≈ -0.1244, M3 ≈ 0.0
#[allow(dead_code)]
pub fn weight_function_edge_crack(x: f64, a: f64, m1: f64, m2: f64, m3: f64) -> f64 {
    if x >= a {
        return 0.0;
    }
    let eta = 1.0 - x / a;
    let denom = (2.0 * PI * (a - x)).sqrt();
    (2.0 / denom) * (1.0 + m1 * eta.sqrt() + m2 * eta + m3 * eta.powf(1.5))
}
/// Compute SIF using the weight function method with numerical integration.
///
/// K = ∫₀ᵃ σ(x) * m(x, a) dx
///
/// Uses the trapezoidal rule with `n_points` integration points.
#[allow(dead_code)]
pub fn sif_weight_function(
    stress_profile: &dyn Fn(f64) -> f64,
    crack_length: f64,
    m1: f64,
    m2: f64,
    m3: f64,
    n_points: usize,
) -> f64 {
    if crack_length <= 0.0 || n_points < 2 {
        return 0.0;
    }
    let dx = crack_length / (n_points - 1) as f64;
    let mut k = 0.0;
    for i in 0..n_points {
        let x = i as f64 * dx;
        let x_safe = x.min(crack_length * (1.0 - 1e-10));
        let sigma = stress_profile(x_safe);
        let m = weight_function_edge_crack(x_safe, crack_length, m1, m2, m3);
        let weight = if i == 0 || i == n_points - 1 {
            0.5
        } else {
            1.0
        };
        k += weight * sigma * m * dx;
    }
    k
}
/// Irwin plasticity correction for the crack-tip plastic zone size.
///
/// Plane stress: r_p = (1/(2π)) * (K_I / σ_ys)²
/// Plane strain: r_p = (1/(6π)) * (K_I / σ_ys)²
///
/// Returns the plastic zone radius r_p.
#[allow(dead_code)]
pub fn irwin_plastic_zone(k_i: f64, sigma_ys: f64, plane_strain: bool) -> f64 {
    let factor = if plane_strain {
        1.0 / (6.0 * PI)
    } else {
        1.0 / (2.0 * PI)
    };
    factor * (k_i / sigma_ys).powi(2)
}
/// Irwin effective crack length accounting for plasticity.
///
/// a_eff = a + r_p
///
/// The effective crack length gives a corrected SIF:
/// K_eff = σ * √(π * a_eff)
#[allow(dead_code)]
pub fn irwin_effective_crack_length(
    crack_length: f64,
    k_i: f64,
    sigma_ys: f64,
    plane_strain: bool,
) -> f64 {
    crack_length + irwin_plastic_zone(k_i, sigma_ys, plane_strain)
}
/// Iterative Irwin plasticity correction.
///
/// Iterates: a_eff = a + r_p(K(a_eff)) until convergence.
/// Returns the converged effective SIF.
#[allow(dead_code)]
pub fn irwin_corrected_sif(
    stress: f64,
    crack_length: f64,
    sigma_ys: f64,
    plane_strain: bool,
    max_iter: usize,
) -> f64 {
    let mut a_eff = crack_length;
    for _ in 0..max_iter {
        let k = stress * (PI * a_eff).sqrt();
        let rp = irwin_plastic_zone(k, sigma_ys, plane_strain);
        let a_new = crack_length + rp;
        if (a_new - a_eff).abs() < 1e-12 * crack_length {
            break;
        }
        a_eff = a_new;
    }
    stress * (PI * a_eff).sqrt()
}
/// Dugdale (strip-yield) plastic zone size for a center crack.
///
/// c = a * sec(π*σ/(2*σ_ys)) - a
///
/// where c is the plastic zone extent, a is the half-crack length.
#[allow(dead_code)]
pub fn dugdale_plastic_zone(half_crack: f64, stress: f64, sigma_ys: f64) -> f64 {
    let sec_val = 1.0 / (PI * stress / (2.0 * sigma_ys)).cos();
    half_crack * (sec_val - 1.0)
}
/// Dugdale CTOD (crack-tip opening displacement) for a center crack.
///
/// δ = (8 * σ_ys * a) / (π * E) * ln(sec(π*σ/(2*σ_ys)))
#[allow(dead_code)]
pub fn dugdale_ctod(half_crack: f64, stress: f64, sigma_ys: f64, e_modulus: f64) -> f64 {
    let sec_val = 1.0 / (PI * stress / (2.0 * sigma_ys)).cos();
    (8.0 * sigma_ys * half_crack) / (PI * e_modulus) * sec_val.ln()
}
/// CTOD (Crack Tip Opening Displacement) from K and material properties.
///
/// Plane stress: δ = K² / (E * σ_ys)
/// Plane strain: δ = K² / (E * σ_ys) * (1 - ν²)
///
/// The factor m ≈ 1 for plane stress, m ≈ 1-ν² for plane strain is absorbed.
#[allow(dead_code)]
pub fn ctod_from_k(k_i: f64, e_modulus: f64, sigma_ys: f64, nu: f64, plane_strain: bool) -> f64 {
    let factor = if plane_strain { 1.0 - nu * nu } else { 1.0 };
    factor * k_i.powi(2) / (e_modulus * sigma_ys)
}
/// Check CTOD fracture criterion.
///
/// Fracture occurs when CTOD >= CTOD_critical.
#[allow(dead_code)]
pub fn ctod_criterion(ctod: f64, ctod_critical: f64) -> bool {
    ctod >= ctod_critical
}
/// CTOA (Crack Tip Opening Angle) from CTOD and measurement distance.
///
/// ψ = 2 * arctan(δ / (2 * d))
///
/// where δ is the CTOD and d is the distance behind the crack tip
/// at which the opening is measured.
#[allow(dead_code)]
pub fn ctoa_from_ctod(ctod: f64, distance_behind_tip: f64) -> f64 {
    2.0 * (ctod / (2.0 * distance_behind_tip)).atan()
}
/// Check CTOA fracture criterion.
///
/// Stable tearing occurs while CTOA > CTOA_critical.
/// When CTOA drops below the critical value, fracture instability is reached.
#[allow(dead_code)]
pub fn ctoa_criterion(ctoa: f64, ctoa_critical: f64) -> bool {
    ctoa >= ctoa_critical
}
/// GTN yield function Φ.
///
/// Φ = (σ_eq/σ_y)² + 2·q1·f*·cosh(3·q2·σ_m / (2·σ_y)) − (1 + q3·f*²)
///
/// Φ < 0: elastic; Φ = 0: yielding.
#[allow(dead_code)]
pub fn gtn_yield_function(
    sigma_eq: f64,
    sigma_m: f64,
    sigma_y: f64,
    f_star: f64,
    params: &GtnParams,
) -> f64 {
    let term1 = (sigma_eq / sigma_y).powi(2);
    let term2 = 2.0 * params.q1 * f_star * (3.0 * params.q2 * sigma_m / (2.0 * sigma_y)).cosh();
    let term3 = 1.0 + params.q3 * f_star.powi(2);
    term1 + term2 - term3
}
/// Void nucleation rate via a Gaussian plastic-strain-controlled model.
///
/// df_n = (fn_void / (sn·√(2π))) · exp(−½·((ε_p − ε_n)/s_n)²) · dε_p
#[allow(dead_code)]
pub fn nucleation_rate(eps_p: f64, deps_p: f64, fn_void: f64, en: f64, sn: f64) -> f64 {
    let exponent = -0.5 * ((eps_p - en) / sn).powi(2);
    (fn_void / (sn * (2.0 * PI).sqrt())) * exponent.exp() * deps_p
}
/// Void growth rate from plastic volumetric strain increment.
///
/// df_growth = (1 − f) · dε_vol_p
#[allow(dead_code)]
pub fn void_growth_rate(f: f64, deps_vol_p: f64) -> f64 {
    (1.0 - f) * deps_vol_p
}
/// Effective (accelerated) void fraction f* accounting for coalescence.
///
/// f*(f) = f                               if f ≤ fc
///        = fc + (1/q1 − fc)/(ff−fc)*(f−fc) if f > fc
#[allow(dead_code)]
pub fn effective_void_fraction(f: f64, params: &GtnParams) -> f64 {
    if f <= params.fc {
        f
    } else {
        let delta_f_star = (1.0 / params.q1 - params.fc) / (params.ff - params.fc);
        params.fc + delta_f_star * (f - params.fc)
    }
}
/// Perform a single GTN time step: update void fraction.
///
/// Combines growth and nucleation contributions.
#[allow(dead_code)]
pub fn gtn_step(f: f64, eps_p: f64, deps_p: f64, deps_vol_p: f64, params: &GtnParams) -> f64 {
    let df_growth = void_growth_rate(f, deps_vol_p);
    let df_nucleation = nucleation_rate(eps_p, deps_p, params.fn_void, params.en, params.sn);
    (f + df_growth + df_nucleation).clamp(0.0, 1.0)
}
/// Convert between K_I and energy release rate G.
///
/// Plane stress: G = K² / E
/// Plane strain: G = K² * (1 - ν²) / E
#[allow(dead_code)]
pub fn k_to_g(k_i: f64, e_modulus: f64, nu: f64, plane_strain: bool) -> f64 {
    let factor = if plane_strain { 1.0 - nu * nu } else { 1.0 };
    factor * k_i.powi(2) / e_modulus
}
/// Convert from energy release rate G to K_I.
///
/// K_I = √(G * E / factor)
#[allow(dead_code)]
pub fn g_to_k(g: f64, e_modulus: f64, nu: f64, plane_strain: bool) -> f64 {
    let factor = if plane_strain { 1.0 - nu * nu } else { 1.0 };
    (g * e_modulus / factor).sqrt()
}
/// Paris law crack growth rate.
///
/// da/dN = C * (ΔK)^m
///
/// where C and m are material constants.
#[allow(dead_code)]
pub fn paris_law_growth_rate(delta_k: f64, c: f64, m: f64) -> f64 {
    if delta_k <= 0.0 {
        return 0.0;
    }
    c * delta_k.powf(m)
}
/// Paris law with threshold and fracture toughness cutoffs.
///
/// da/dN = C * (ΔK)^m * (1 - ΔK_th/ΔK)^p / (1 - K_max/K_Ic)^q
///
/// (Forman/Mettu NASGRO equation simplified form)
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn paris_forman_growth_rate(
    delta_k: f64,
    k_max: f64,
    c: f64,
    m: f64,
    delta_k_th: f64,
    k_ic: f64,
) -> f64 {
    if delta_k <= delta_k_th || k_max >= k_ic {
        return if k_max >= k_ic { f64::INFINITY } else { 0.0 };
    }
    let threshold_factor = 1.0 - delta_k_th / delta_k;
    let instability_factor = 1.0 - k_max / k_ic;
    if instability_factor <= 0.0 {
        return f64::INFINITY;
    }
    c * delta_k.powf(m) * threshold_factor / instability_factor
}
/// Mode I SIF for a single-edge-notched beam (SENB) in three-point bending.
///
/// K_I = F * S / (B * W^0.5) * f(a/W)
///
/// where F is the applied force, S is the support span, B is the thickness,
/// W is the beam depth (height), a is the crack length, and the geometry
/// function f(α) follows the ASTM E1820 / Srawley polynomial:
///
/// f(α) = (3*α^0.5 * (1.99 - α*(1-α)*(2.15 - 3.93*α + 2.7*α²)))
///         / (2*(1+2*α)*(1-α)^1.5)
///
/// Valid for α = a/W in \[0.0, 1.0\].
#[allow(dead_code)]
pub fn sif_three_point_bend(
    force: f64,
    span: f64,
    thickness: f64,
    width: f64,
    crack_length: f64,
) -> f64 {
    let alpha = crack_length / width;
    let alpha_sqrt = alpha.sqrt();
    let poly = 1.99 - alpha * (1.0 - alpha) * (2.15 - 3.93 * alpha + 2.7 * alpha * alpha);
    let f_alpha = (3.0 * alpha_sqrt * poly) / (2.0 * (1.0 + 2.0 * alpha) * (1.0 - alpha).powf(1.5));
    (force * span) / (thickness * width.sqrt()) * f_alpha
}
/// Critical crack length for a SENB specimen at fracture.
///
/// Solves K_I(a) = K_Ic for a given stress state using bisection.
/// Returns None if no valid root is found in \[a_min, a_max\].
#[allow(dead_code)]
pub fn senb_critical_crack_length(
    force: f64,
    span: f64,
    thickness: f64,
    width: f64,
    k_ic: f64,
    n_iter: usize,
) -> Option<f64> {
    let a_min = 0.01 * width;
    let a_max = 0.9 * width;
    let f_lo = sif_three_point_bend(force, span, thickness, width, a_min) - k_ic;
    let f_hi = sif_three_point_bend(force, span, thickness, width, a_max) - k_ic;
    if f_lo * f_hi > 0.0 {
        return None;
    }
    let mut lo = a_min;
    let mut hi = a_max;
    for _ in 0..n_iter {
        let mid = 0.5 * (lo + hi);
        let f_mid = sif_three_point_bend(force, span, thickness, width, mid) - k_ic;
        if f_mid * f_lo < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
        if (hi - lo) / width < 1e-12 {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}
/// Mixed-mode SIF for an inclined through-crack in an infinite plate.
///
/// For a crack inclined at angle β to the loading direction:
/// K_I  = σ * √(π*a) * cos²(β)
/// K_II = σ * √(π*a) * sin(β) * cos(β)
///
/// Returns (K_I, K_II).
#[allow(dead_code)]
pub fn sif_mixed_mode_inclined_crack(stress: f64, half_crack: f64, angle_rad: f64) -> (f64, f64) {
    let k_base = stress * (PI * half_crack).sqrt();
    let k_i = k_base * angle_rad.cos().powi(2);
    let k_ii = k_base * angle_rad.sin() * angle_rad.cos();
    (k_i, k_ii)
}
/// Mode III (anti-plane shear) SIF for an edge crack.
///
/// K_III = τ * √(π * a) where τ is the shear stress.
#[allow(dead_code)]
pub fn sif_mode_iii_edge_crack(shear_stress: f64, crack_length: f64) -> f64 {
    shear_stress * (PI * crack_length).sqrt()
}
/// Effective (equivalent) SIF for mixed-mode I+II loading.
///
/// K_eff = sqrt(K_I² + K_II²)
///
/// A common approximation for the effective driving force for crack propagation
/// in mixed-mode loading.
#[allow(dead_code)]
pub fn sif_effective_mixed_mode(k_i: f64, k_ii: f64) -> f64 {
    (k_i.powi(2) + k_ii.powi(2)).sqrt()
}
/// Effective SIF including mode III under plane strain.
///
/// K_eff = sqrt(K_I² + K_II² + K_III²/(1-ν))
#[allow(dead_code)]
pub fn sif_effective_mixed_mode_iii(k_i: f64, k_ii: f64, k_iii: f64, nu: f64) -> f64 {
    (k_i.powi(2) + k_ii.powi(2) + k_iii.powi(2) / (1.0 - nu)).sqrt()
}
/// Maximum hoop stress angle for mixed mode I+II.
///
/// The crack propagates in the direction θ_c that maximises the hoop stress:
/// θ_c = 2 * arctan(K_I/(4*K_II) - sqrt((K_I/(4*K_II))² + 0.5))
///
/// Returns angle in radians (θ = 0 for pure mode I).
#[allow(dead_code)]
pub fn mixed_mode_propagation_angle(k_i: f64, k_ii: f64) -> f64 {
    if k_ii.abs() < 1e-30 {
        return 0.0;
    }
    let ratio = k_i / (4.0 * k_ii);
    let inner = (ratio * ratio + 0.5).sqrt();
    2.0 * (ratio - inner).atan()
}
/// Walker equation for R-ratio correction of fatigue crack growth.
///
/// Computes the effective ΔK accounting for the stress ratio R:
///
/// ΔK_eff = ΔK / (1 - R)^(1 - γ)
///
/// where γ is the Walker exponent (0 < γ < 1, typically 0.5 for metals).
/// At R = 0 and γ = 0.5 this reduces to ΔK_eff = ΔK / √(1-R) = √(ΔK * K_max).
#[allow(dead_code)]
pub fn walker_delta_k_eff(delta_k: f64, r_ratio: f64, gamma: f64) -> f64 {
    let denom = (1.0 - r_ratio).powf(1.0 - gamma);
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    delta_k / denom
}
/// Paris law crack growth rate with Walker R-ratio correction.
///
/// da/dN = C * (ΔK_eff)^m = C * (ΔK / (1-R)^(1-γ))^m
#[allow(dead_code)]
pub fn paris_walker_rate(delta_k: f64, r_ratio: f64, c: f64, m: f64, gamma: f64) -> f64 {
    if delta_k <= 0.0 {
        return 0.0;
    }
    let dk_eff = walker_delta_k_eff(delta_k, r_ratio, gamma);
    c * dk_eff.powf(m)
}
/// Forman equation for crack growth with R-ratio.
///
/// da/dN = C * (ΔK)^m / ((1-R) * K_c - ΔK)
///
/// Accounts for both mean stress ratio effects and fracture instability.
#[allow(dead_code)]
pub fn forman_crack_growth_rate(delta_k: f64, r_ratio: f64, c: f64, m: f64, k_c: f64) -> f64 {
    if delta_k <= 0.0 {
        return 0.0;
    }
    let denom = (1.0 - r_ratio) * k_c - delta_k;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    c * delta_k.powf(m) / denom
}
/// First-order Irwin plastic zone size.
///
/// Plane stress: r_p = K_I² / (2π σ_ys²)
/// Plane strain: r_p = K_I² / (6π σ_ys²)
///
/// (Same as the existing `irwin_plastic_zone` — provided here for symmetry with
/// the second-order version below and to allow standalone use.)
#[allow(dead_code)]
pub fn plastic_zone_size_first_order(k_i: f64, sigma_ys: f64, plane_strain: bool) -> f64 {
    let factor = if plane_strain {
        1.0 / (6.0 * PI)
    } else {
        1.0 / (2.0 * PI)
    };
    factor * (k_i / sigma_ys).powi(2)
}
/// Second-order (iterative) Irwin plastic zone size including plasticity correction.
///
/// Iterates: r_p^(n+1) = (1/(2π)) * (K_I(a + r_p^(n)) / σ_ys)²
///
/// until convergence, accounting for the increase in SIF due to the plastic zone.
#[allow(dead_code)]
pub fn plastic_zone_size_iterative(
    stress: f64,
    crack_length: f64,
    sigma_ys: f64,
    plane_strain: bool,
    max_iter: usize,
) -> f64 {
    let factor = if plane_strain {
        1.0 / (6.0 * PI)
    } else {
        1.0 / (2.0 * PI)
    };
    let mut rp = 0.0_f64;
    for _ in 0..max_iter {
        let a_eff = crack_length + rp;
        let k = stress * (PI * a_eff).sqrt();
        let rp_new = factor * (k / sigma_ys).powi(2);
        if (rp_new - rp).abs() < 1e-15 * (1.0 + rp) {
            return rp_new;
        }
        rp = rp_new;
    }
    rp
}
/// Constraint factor for the plastic zone accounting for biaxiality.
///
/// Under biaxial loading with biaxiality ratio λ = σ_2 / σ_1,
/// the effective yield stress is modified:
///
/// σ_ys_eff = σ_ys / √(1 - λ + λ²)  (von Mises)
///
/// The biaxiality-corrected plastic zone is r_p = K_I²/(2π σ_ys_eff²).
#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn plastic_zone_biaxial(k_i: f64, sigma_ys: f64, lambda: f64, plane_strain: bool) -> f64 {
    let vm_factor = (1.0 - lambda + lambda * lambda).sqrt();
    let sigma_ys_eff = sigma_ys / vm_factor;
    let factor = if plane_strain {
        1.0 / (6.0 * PI)
    } else {
        1.0 / (2.0 * PI)
    };
    factor * (k_i / sigma_ys_eff).powi(2)
}
/// Irwin CTOD formula: δ = K_I² / (σ_ys * E')
///
/// where E' = E (plane stress) or E' = E / (1 - ν²) (plane strain).
#[allow(dead_code)]
pub fn ctod_irwin(k_i: f64, sigma_ys: f64, e: f64, nu: f64, plane_strain: bool) -> f64 {
    let e_prime = if plane_strain { e / (1.0 - nu * nu) } else { e };
    k_i.powi(2) / (sigma_ys * e_prime)
}
/// Wells CTOD formula with Irwin plastic zone correction.
///
/// δ = (4/π) * (K_I² / (E * σ_ys))  (plane stress approximation)
#[allow(dead_code)]
pub fn ctod_wells(k_i: f64, e: f64, sigma_ys: f64) -> f64 {
    (4.0 / PI) * k_i.powi(2) / (e * sigma_ys)
}
/// CTOD from J-integral (equivalence relation).
///
/// δ = J / (m * σ_ys)
///
/// where m ≈ 1 (plane stress) or m ≈ 2 (plane strain).
#[allow(dead_code)]
pub fn ctod_from_j(j_integral: f64, sigma_ys: f64, plane_strain: bool) -> f64 {
    let m = if plane_strain { 2.0 } else { 1.0 };
    j_integral / (m * sigma_ys)
}
/// J-integral from CTOD.
///
/// J = m * σ_ys * δ
#[allow(dead_code)]
pub fn j_from_ctod(ctod: f64, sigma_ys: f64, plane_strain: bool) -> f64 {
    let m = if plane_strain { 2.0 } else { 1.0 };
    m * sigma_ys * ctod
}
/// Minimum specimen thickness for valid plane-strain fracture toughness test.
///
/// Per ASTM E399: B_min = 2.5 * (K_Ic / σ_ys)²
#[allow(dead_code)]
pub fn minimum_thickness_plane_strain(k_ic: f64, sigma_ys: f64) -> f64 {
    2.5 * (k_ic / sigma_ys).powi(2)
}
/// Check whether the test satisfies plane-strain validity per ASTM E399.
///
/// Valid if: B ≥ B_min AND a ≥ B_min AND (W - a) ≥ B_min
#[allow(dead_code)]
pub fn is_plane_strain_valid(
    thickness: f64,
    crack_length: f64,
    width: f64,
    k_ic: f64,
    sigma_ys: f64,
) -> bool {
    let b_min = minimum_thickness_plane_strain(k_ic, sigma_ys);
    thickness >= b_min && crack_length >= b_min && (width - crack_length) >= b_min
}
/// Convert K_Ic to critical J (Jc) using the plane-strain formula.
///
/// J_c = K_Ic² * (1 - ν²) / E
#[allow(dead_code)]
pub fn k_to_j(k_ic: f64, e: f64, nu: f64) -> f64 {
    k_ic.powi(2) * (1.0 - nu * nu) / e
}
/// Convert critical J to K_Ic.
///
/// K_Ic = sqrt(J_c * E / (1 - ν²))
#[allow(dead_code)]
pub fn j_to_k(j_c: f64, e: f64, nu: f64) -> f64 {
    (j_c * e / (1.0 - nu * nu)).sqrt()
}
/// Dynamic stress intensity factor (simple Freund approximation).
///
/// K_dyn = K_static * (1 - v/c_R)^(1/2)
///
/// where v is the crack velocity and c_R is the Rayleigh wave speed.
/// At v = 0, K_dyn = K_static. At v = c_R, K_dyn = 0 (terminal velocity).
#[allow(dead_code)]
pub fn sif_dynamic_freund(k_static: f64, crack_velocity: f64, rayleigh_speed: f64) -> f64 {
    if rayleigh_speed <= 0.0 || crack_velocity >= rayleigh_speed {
        return 0.0;
    }
    let ratio = crack_velocity / rayleigh_speed;
    k_static * (1.0 - ratio).sqrt()
}
/// Rayleigh wave speed (approximate).
///
/// c_R ≈ c_s * (0.862 + 1.14*ν) / (1 + ν)
///
/// where c_s = sqrt(G/ρ) is the shear wave speed.
#[allow(dead_code)]
pub fn rayleigh_wave_speed(shear_modulus: f64, density: f64, nu: f64) -> f64 {
    let c_s = (shear_modulus / density).sqrt();
    c_s * (0.862 + 1.14 * nu) / (1.0 + nu)
}
/// Crack arrest condition: K_static ≤ K_Ia (arrest toughness).
///
/// Returns true if the crack will arrest (K_static < K_Ia).
#[allow(dead_code)]
pub fn will_crack_arrest(k_static: f64, k_ia: f64) -> bool {
    k_static < k_ia
}
/// Mode mixity ratio ψ = arctan(K_II / K_I) in radians.
///
/// ψ = 0 → pure mode I; ψ = π/2 → pure mode II.
#[allow(dead_code)]
pub fn mode_mixity_angle(k_i: f64, k_ii: f64) -> f64 {
    k_ii.atan2(k_i)
}
/// Mode mixity phase angle (degrees).
#[allow(dead_code)]
pub fn mode_mixity_angle_degrees(k_i: f64, k_ii: f64) -> f64 {
    mode_mixity_angle(k_i, k_ii).to_degrees()
}
/// Fracture mode I fraction from K_I, K_II, K_III.
///
/// η_I = K_I² / (K_I² + K_II² + K_III²)
#[allow(dead_code)]
pub fn mode_i_fraction(k_i: f64, k_ii: f64, k_iii: f64) -> f64 {
    let total = k_i * k_i + k_ii * k_ii + k_iii * k_iii;
    if total <= 0.0 {
        return 1.0;
    }
    k_i * k_i / total
}
/// Stress triaxiality T = σ_m / σ_eq.
///
/// Used in ductile fracture models (GTN, Rice-Tracey, etc.).
/// High triaxiality accelerates void growth and ductile fracture.
#[allow(dead_code)]
pub fn stress_triaxiality(sigma_m: f64, sigma_eq: f64) -> f64 {
    if sigma_eq.abs() < 1e-30 {
        return 0.0;
    }
    sigma_m / sigma_eq
}
/// Rice-Tracey void growth factor.
///
/// R_dot/R = C * exp(3*σ_m / (2*σ_0)) * eps_dot_eq
///
/// C ≈ 0.283 (original R-T value).
#[allow(dead_code)]
pub fn rice_tracey_void_growth(sigma_m: f64, sigma_0: f64, eps_dot_eq: f64) -> f64 {
    let c = 0.283;
    c * (1.5 * sigma_m / sigma_0).exp() * eps_dot_eq
}
/// Stress constraint factor T* (T-stress normalised by yield stress).
///
/// T* = T_stress / σ_ys
///
/// Accounts for in-plane constraint effects on fracture toughness.
#[allow(dead_code)]
pub fn t_stress_constraint(t_stress: f64, sigma_ys: f64) -> f64 {
    if sigma_ys.abs() < 1e-30 {
        return 0.0;
    }
    t_stress / sigma_ys
}
/// Coffin-Manson low-cycle fatigue (LCF) relation.
///
/// Δεp / 2 = εf' * (2*Nf)^c
///
/// where Δεp is the plastic strain range, εf' is the fatigue ductility
/// coefficient, and c is the fatigue ductility exponent (typically -0.5 to -0.7).
#[allow(dead_code)]
pub fn coffin_manson_strain_amplitude(n_f: f64, eps_f_prime: f64, c: f64) -> f64 {
    eps_f_prime * (2.0 * n_f).powf(c)
}
/// Basquin high-cycle fatigue (HCF) relation.
///
/// Δσ/2 = σf' * (2*Nf)^b
///
/// where σf' is the fatigue strength coefficient, b is the fatigue strength
/// exponent (typically -0.05 to -0.12).
#[allow(dead_code)]
pub fn basquin_stress_amplitude(n_f: f64, sigma_f_prime: f64, b: f64) -> f64 {
    sigma_f_prime * (2.0 * n_f).powf(b)
}
/// Total strain amplitude from combined Coffin-Manson + Basquin (Morrow).
///
/// Δε/2 = σf'/E * (2*Nf)^b + εf' * (2*Nf)^c
#[allow(dead_code)]
pub fn morrow_total_strain_amplitude(
    n_f: f64,
    sigma_f_prime: f64,
    b: f64,
    eps_f_prime: f64,
    c: f64,
    e_modulus: f64,
) -> f64 {
    sigma_f_prime / e_modulus * (2.0 * n_f).powf(b) + eps_f_prime * (2.0 * n_f).powf(c)
}
/// Fatigue life (cycles to failure) from stress amplitude using Basquin.
///
/// Nf = (1/2) * (Δσ/(2*σf'))^(1/b)
#[allow(dead_code)]
pub fn basquin_fatigue_life(delta_sigma_over_2: f64, sigma_f_prime: f64, b: f64) -> f64 {
    if sigma_f_prime <= 0.0 || b >= 0.0 {
        return f64::INFINITY;
    }
    0.5 * (delta_sigma_over_2 / sigma_f_prime).powf(1.0 / b)
}
/// Smith-Watson-Topper (SWT) parameter for mean stress correction.
///
/// P_SWT = σ_max * Δε/2
///
/// For mean-stress corrected fatigue life.
#[allow(dead_code)]
pub fn swt_parameter(sigma_max: f64, strain_amplitude: f64) -> f64 {
    sigma_max * strain_amplitude
}
/// Mode I SIF for a through crack emanating from a circular hole in a plate.
///
/// K_I = σ * √(π*a) * F(λ)
///
/// where λ = a/R (crack length to hole radius ratio), and F(λ) ≈ 1 + 0.1215*λ − 0.0598*λ² (approximate).
#[allow(dead_code)]
pub fn sif_hole_edge_crack(stress: f64, crack_length: f64, hole_radius: f64) -> f64 {
    let lambda = crack_length / hole_radius;
    let f = 1.0 + 0.1215 * lambda - 0.0598 * lambda * lambda;
    stress * (PI * crack_length).sqrt() * f.max(1.0)
}
/// Mode I SIF for a surface semi-elliptical crack in a plate.
///
/// K_I = σ * √(π*a/Q) * F_s
///
/// where Q is the shape factor (Q ≈ φ² - 0.212*(σ/σ_ys)², φ = elliptic integral).
/// Simplified: Q ≈ (1 + 4.593*(a/c)^1.65) for a/c ≤ 1.
#[allow(dead_code)]
pub fn sif_surface_crack(stress: f64, depth: f64, half_length: f64) -> f64 {
    let ratio = depth / half_length;
    let phi2 = if ratio <= 1.0 {
        1.0 + 4.593 * ratio.powf(1.65)
    } else {
        1.0 + 4.593 * (1.0 / ratio).powf(1.65)
    };
    stress * (PI * depth / phi2).sqrt()
}
/// Near-tip stress field (mode I) in polar coordinates.
///
/// σ_rr = K_I / √(2π*r) * cos(θ/2) * \[1 - sin(θ/2)*sin(3θ/2)\]
/// σ_θθ = K_I / √(2π*r) * cos(θ/2) * \[1 + sin(θ/2)*sin(3θ/2)\]
/// σ_rθ = K_I / √(2π*r) * sin(θ/2) * cos(θ/2) * cos(3θ/2)
///
/// Returns (sigma_rr, sigma_theta_theta, sigma_r_theta).
#[allow(dead_code)]
pub fn near_tip_stress_mode_i(k_i: f64, r: f64, theta: f64) -> (f64, f64, f64) {
    if r <= 0.0 {
        return (f64::INFINITY, f64::INFINITY, f64::INFINITY);
    }
    let coeff = k_i / (2.0 * PI * r).sqrt();
    let half = theta / 2.0;
    let sigma_rr = coeff * half.cos() * (1.0 - half.sin() * (3.0 * half).sin());
    let sigma_tt = coeff * half.cos() * (1.0 + half.sin() * (3.0 * half).sin());
    let sigma_rt = coeff * half.sin() * half.cos() * (3.0 * half).cos();
    (sigma_rr, sigma_tt, sigma_rt)
}
/// Numerical integration of Paris law da/dN = C*(ΔK)^m.
///
/// Integrates from initial crack length a0 to final af (or until K_max ≥ K_Ic).
/// Returns total number of cycles.
#[allow(dead_code)]
pub fn paris_law_cycles(
    a0: f64,
    a_f: f64,
    stress_range: f64,
    c: f64,
    m: f64,
    k_ic: f64,
    n_steps: usize,
) -> f64 {
    let mut a = a0;
    let mut n_total = 0.0;
    let da = (a_f - a0) / n_steps as f64;
    for _ in 0..n_steps {
        let k_max = stress_range * (PI * a).sqrt();
        if k_max >= k_ic {
            break;
        }
        let delta_k = stress_range * (PI * a).sqrt();
        let da_dn = paris_law_growth_rate(delta_k, c, m);
        if da_dn <= 0.0 {
            break;
        }
        n_total += da / da_dn;
        a += da;
        if a >= a_f {
            break;
        }
    }
    n_total
}
/// Forman equation fatigue life integration.
///
/// Integrates da/dN = C*(ΔK)^m / ((1-R)*K_c - ΔK).
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn forman_law_cycles(
    a0: f64,
    a_f: f64,
    delta_sigma: f64,
    r_ratio: f64,
    c: f64,
    m: f64,
    k_c: f64,
    n_steps: usize,
) -> f64 {
    let mut a = a0;
    let mut n_total = 0.0;
    let da = (a_f - a0) / n_steps as f64;
    for _ in 0..n_steps {
        let delta_k = delta_sigma * (PI * a).sqrt();
        let k_max = delta_k / (1.0 - r_ratio).max(1e-10);
        if k_max >= k_c {
            break;
        }
        let da_dn = forman_crack_growth_rate(delta_k, r_ratio, c, m, k_c);
        if da_dn <= 0.0 || !da_dn.is_finite() {
            break;
        }
        n_total += da / da_dn;
        a += da;
        if a >= a_f {
            break;
        }
    }
    n_total
}
